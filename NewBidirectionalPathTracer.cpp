#include "StdAfx.h"
#include "NewBidirectionalPathTracer.h"

static FILE* fp = fopen("debug_new_bpt.txt" , "w");

vector<vec3f> NewBidirectionalPathTracer::renderPixels(const Camera& camera)
{
	unsigned t_start = clock();

	vector<vec3f> pixelColors(camera.width * camera.height, vec3f(0, 0, 0));
	vector<omp_lock_t> pixelLocks(pixelColors.size());

	preprocessEmissionSampler();

	for(int i=0; i<pixelLocks.size(); i++)
	{
		omp_init_lock(&pixelLocks[i]);
	}

	omp_lock_t cmdLock;
	omp_init_lock(&cmdLock);

	for(unsigned s=0; s<spp; s++)
	{
		vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));

		string cmd;

		unsigned t = clock();

		if(!renderer->scene.usingGPU())
		{
			lightStateIndex.resize(lightPathNum);
			memset(&lightStateIndex[0] , 0 , lightStateIndex.size() * sizeof(int));

			lightStates.reserve(lightPathNum);
			lightStates.clear();
			cameraStates.reserve(cameraPathNum);
			cameraStates.clear();

			for(int p=0; p<lightPathNum; p++)
			{
				Path lightPath;

				Ray lightRay = genEmissiveSurfaceSample(true , false);

				samplePath(lightPath, lightRay);

				BidirPathState lightState;
				int s = 1;

				genLightSample(lightPath , lightState);

				//fprintf(fp , "==================\n");
				for (s = 1; s < lightPath.size(); s++)
				{
					lightState.pos = lightPath[s].origin;
					lightState.dir = lightPath[s].direction;
					float dist = (lightPath[s].origin - lightPath[s - 1].origin).length();
					if (abs(dist) < 1e-6f)
						break;
					vec3f decayFactor = lightPath[s - 1].getRadianceDecay(dist);
					lightState.throughput *= decayFactor;

					//fprintf(fp , "thr=(%.8f,%.8f,%.8f) , cos=%.8f, pdf=%.8f\n" , 
					//	lightState.throughput.x , lightState.throughput.y , lightState.throughput.z ,
					//	lightPath[s].getCosineTerm() , lightPath[s].directionProb);

					if (lightState.pathLength > 1 || lightState.isFiniteLight)
						lightState.dVCM *= mis(dist * dist);
					float cosWi = lightPath[s].getContactNormal().dot(-lightPath[s].direction);
					lightState.dVCM /= mis(abs(cosWi));
					lightState.dVC /= mis(abs(cosWi));

					if (lightPath[s].directionSampleType == Ray::RANDOM)
						lightStates.push_back(lightState);

					if (lightPath[s].directionSampleType == Ray::RANDOM)
					{
						int _x(0) , _y(0);
						vec3f color(0.f);
						color = colorByConnectingCamera(camera , lightState , lightPath[s] , lightPath[s - 1] , _x , _y);
				
						if (y(color) > 0)
						{
							omp_set_lock(&pixelLocks[_y*camera.width + _x]);
							singleImageColors[_y*camera.width + _x] += color;
							omp_unset_lock(&pixelLocks[_y*camera.width + _x]);
						}
					}

					if (s + 1 >= lightPath.size())
						break;
					
					float bsdfDirPdf = lightPath[s].directionProb;
					Ray inRay = lightPath[s + 1];
					inRay.direction = -lightPath[s].direction;
					Ray outRay = lightPath[s];
					outRay.direction = -lightPath[s - 1].direction;
					float bsdfRevPdf = inRay.getDirectionSampleProbDensity(outRay);
					
					lightState.throughput = (lightState.throughput * lightPath[s].color) *
						(lightPath[s].getCosineTerm() / lightPath[s].directionProb);

					if (lightPath[s].directionSampleType == Ray::DEFINITE)
					{
						lightState.dVCM = 0.f;
						if (abs(bsdfRevPdf - lightPath[s].directionProb) > 1e-6f) 
						{
							printf("error: dir = %.8f, rev = %.8f\n" , lightPath[s].directionProb , 
								bsdfRevPdf);
							bsdfRevPdf = lightPath[s].directionProb;
						}
						lightState.dVC *= mis(lightPath[s].getCosineTerm());
						lightState.specularPath &= 1;
					}
					else
					{
						lightState.dVC = mis(lightPath[s].getCosineTerm() / lightPath[s].directionProb) *
							(lightState.dVC * mis(bsdfRevPdf) + lightState.dVCM);
						lightState.dVCM = mis(1.f / lightPath[s].directionProb);
						lightState.specularPath &= 0;
					}
				}
				lightStateIndex[p] = (int)lightStates.size();
			}
			cameraPathNum = 0;
			for (int p=0; p<cameraPathNum; p++)
			{
				Path cameraPath;

				Ray cameraRay = camera.generateRay(p);

				samplePath(cameraPath, cameraRay);

				BidirPathState cameraState;
				vec3f color(0.f);

				genCameraSample(camera , cameraPath , cameraState);

				for (int s = 1; s < cameraPath.size(); s++)
				{
					float dist = (cameraPath[s].origin - cameraPath[s - 1].origin).length();
					if (dist < 1e-6f)
						break;
					cameraState.dVCM *= mis(dist * dist);
					float cosWi = cameraPath[s].getContactNormal().dot(-cameraPath[s - 1].direction);
					cameraState.dVCM /= mis(abs(cosWi));
					cameraState.dVC /= mis(abs(cosWi));

					if (cameraPath[s].contactObject && cameraPath[s].contactObject->emissive())
					{

					}
				}
			}

			if(cmd == "exit")
				return pixelColors;

			//eliminateVignetting(singleImageColors);

			for(int i=0; i<pixelColors.size(); i++)
			{
				pixelColors[i] *= s / float(s + 1);
				pixelColors[i] += singleImageColors[i] / (s + 1);//*camera.width*camera.height;
			}

			printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, (clock()-t_start)/1000);

			if (clock() / 1000 >= lastTime)
			{
				showCurrentResult(pixelColors , &lastTime);
				lastTime += timeInterval;
			}
			else
				showCurrentResult(pixelColors);
		}
	}
	return pixelColors;
}

void NewBidirectionalPathTracer::genLightSample(Path& lightPath , BidirPathState& lightState)
{
	lightState.pos = lightPath[0].origin;
	lightState.dir = lightPath[0].direction;
	float emitPdf = lightPath[0].directionProb * lightPath[0].originProb;
	lightState.throughput = lightPath[0].color * lightPath[0].getCosineTerm() / emitPdf;
	lightState.pathLength = 1;
	lightState.isFiniteLight = 1;
	lightState.specularPath = 1;
	lightState.specularVertexNum = 0;
	float dirPdf = lightPath[0].originProb;

	lightState.dVCM = mis(dirPdf / emitPdf);
	lightState.dVC = mis(1.f / emitPdf);
}

vec3f NewBidirectionalPathTracer::colorByConnectingCamera(const Camera& camera, const BidirPathState& lightState , const Ray& ray , const Ray& lastRay , int& _x , int& _y)
{
	vec2<float> pCoord = camera.transToPixel(lightState.pos);
	int x = pCoord.x;
	int y = pCoord.y;
	
	if(!(x >= 0 && x < camera.width && y >= 0 && y < camera.height))
		return vec3f(0.f);
	_x = x; _y = y;

	vec3f dirToCamera = camera.position - lightState.pos;
	vec3f forward = camera.focus - camera.position;
	forward.normalize();

	if (forward.dot(-dirToCamera) <= 0)
		return vec3f(0.f);

	float dist = dirToCamera.length();
	float distEye2 = dist * dist;
	float cameraDistToScreen2 = camera.sightDist * camera.sightDist;
	dirToCamera = dirToCamera / dist;

	float cosToCamera;
	cosToCamera = std::abs(ray.getContactNormal().dot(dirToCamera));
	if (cosToCamera < 1e-6f)
		return vec3f(0.f);
	//printf("cosToCamera = %.8f\n" , cosToCamera);

	Ray outRay = ray;
	outRay.direction = dirToCamera;

	vec3f bsdfFactor = lastRay.getBSDF(outRay);
	if (intensity(bsdfFactor) < 1e-6f)
		return vec3f(0.f);

	float cosAtCamera = forward.dot(-dirToCamera);
	//printf("cosAtCamera = %.8f\n" , cosAtCamera);

	float imagePointToCameraDist = camera.sightDist / cosAtCamera;
	float imageToSolidAngleFactor = imagePointToCameraDist *
		imagePointToCameraDist / cosAtCamera;
	float imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / distEye2;

	float cameraPdfArea = imageToSurfaceFactor * 1.f; // pixel area is 1
	
	float surfaceToImageFactor = 1.f / imageToSurfaceFactor;

	vec3f totContrib = lightState.throughput;

	//---- still buggy, fix me ----
	vec3f color = (totContrib * bsdfFactor) //* cosAtCamera * cosToCamera / distEye2;
		 / (surfaceToImageFactor * lightPathNum);
	//color *= powf(cosAtCamera , 4.f) / pixelNum;
	//-----------------------------

	outRay.direction = -lastRay.direction;
	Ray inRay;
	inRay.origin = camera.position;
	inRay.direction = -dirToCamera;

	float bsdfRevPdf = inRay.getDirectionSampleProbDensity(outRay);

	if (!testVisibility(inRay , outRay))
		return vec3f(0.f);

	float wLight = mis(cameraPdfArea / lightPathNum) * (lightState.dVCM +
		mis(bsdfRevPdf) * lightState.dVC);
	float weight = 1.f / (wLight + 1.f);

	//color *= weight;
	/*
	if (lightState.isSpecularPath && bsdfFactor[2] > bsdfFactor[1] && bsdfFactor[2] > bsdfFactor[0])
	{
		fprintf(fp , "============blue============\n");
		vec3f resx = color * lightPathNum / powf(cosAtCamera , 4.f);
		fprintf(fp , "factor = %.8f, cosToCamera = %.8f, pdf = %.8f , weight = %.6lf,\nbsdfFactor = (%.8f,%.8f,%.8f), connect camera = (%.6lf,%.6lf,%.6lf)\ntotContrib = (%.8f,%.8f,%.8f), thr = (%.8f,%.8f,%.8f)\n" ,
			lightPathNum * surfaceToImageFactor , cosToCamera , pdf , weightFactor , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] ,
			resx[0] , resx[1] , resx[2] , totContrib[0] , totContrib[1] , totContrib[2] , 
			lightState.throughput[0] , lightState.throughput[1] , lightState.throughput[2]);
	}
	*/
	/*
	vec3f resx = color * lightPathNum / powf(cosAtCamera , 4.f);
	fprintf(fp , "cosToCamera = %.8f, pdf = %.8f , weight = %.6lf,\nbsdfFactor = (%.8f,%.8f,%.8f), connect camera = (%.6lf,%.6lf,%.6lf)\n" ,
		cosToCamera , pdf , weightFactor , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] ,
		resx[0] , resx[1] , resx[2]);
	*/

	return color;
}

void NewBidirectionalPathTracer::genCameraSample(const Camera& camera , Path& cameraPath , BidirPathState& cameraState)
{
	cameraState.pos = cameraPath[0].origin;
	cameraState.dir = cameraPath[0].direction;
	cameraState.pathLength = 1;
	cameraState.specularPath = 1;
	cameraState.specularVertexNum = 0;
	// Different!
	cameraState.throughput = vec3f(1.f);
	// =========

	vec3f forward = camera.focus - camera.position;
	forward.normalize();
	float cosAtCamera = cameraPath[0].direction.dot(forward);
	float imagePointToCameraDist = camera.sightDist / cosAtCamera;
	float imageToSolidAngle = imagePointToCameraDist * imagePointToCameraDist / cosAtCamera;

	float cameraPdf = imageToSolidAngle;
	cameraState.dVCM = mis(lightPathNum / cameraPdf);
	cameraState.dVC = 0.f;
}
