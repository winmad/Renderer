#include "StdAfx.h"
#include "IptTracer.h"
#include "SceneEmissiveObject.h"
#include "UniformSphericalSampler.h"
#include "NoSelfIntersectionCondition.h"

static FILE *fp = fopen("debug_ipt.txt" , "w");

float INV_2_PI = 0.5 / M_PI;

vector<vec3f> IptTracer::renderPixels(const Camera& camera)
{
	vector<vec3f> pixelColors(camera.width * camera.height, vec3f(0, 0, 0));
	vector<omp_lock_t> pixelLocks(pixelColors.size());

	preprocessEmissionSampler();
	preprocessOtherSampler();

	totArea = renderer->scene.getTotalArea();

	for(int i=0; i<pixelLocks.size(); i++)
	{
		omp_init_lock(&pixelLocks[i]);
	}

	omp_lock_t cmdLock;
	omp_init_lock(&cmdLock);

	Real r0 = mergeRadius;

	for(unsigned s=0; s<spp; s++)
	{
		if (!renderer->scene.usingGPU())
		{
			partPathMergeIndex.resize(lightPathNum + interPathNum);

			partialSubPathList.clear();

			mergeRadius = r0 * (pow(s+1, 0.5f*(alpha-1)));
			mergeRadius = std::max(mergeRadius , 1e-7f);
			printf("mergeRadius = %.8f\n" , mergeRadius);

			mergeKernel = 1.f / (M_PI * mergeRadius * 
				mergeRadius * (Real)partialPathNum);

			vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));
	
			string cmd;
	
			unsigned t = clock();

			vector<Path*> lightPathList(lightPathNum , NULL);
			vector<Path*> interPathList(interPathNum, NULL);

			genLightPaths(cmdLock , lightPathList);

			genIntermediatePaths(cmdLock , interPathList);
	
			mergePartialPaths(cmdLock);

			printf("partialPathNum = %d\n" , partialSubPathList.size());

#pragma omp parallel for
			for (int i = 0; i < partialSubPathList.size(); i++)
			{
				IptPathState& subPath = partialSubPathList[i];
				colorByConnectingCamera(pixelLocks , camera , singleImageColors , subPath);

				
				if (s == 0)
				{
				fprintf(fp , "dirContrib=(%.4f,%.4f,%.4f), indirContrib=(%.4f,%.4f,%.4f), throughput=(%.4f,%.4f,%.4f)\n" ,
					subPath.dirContrib[0] , subPath.dirContrib[1] , subPath.dirContrib[2] ,
					subPath.indirContrib[0] , subPath.indirContrib[1] , subPath.indirContrib[2] ,
					subPath.throughput[0] , subPath.throughput[1] , subPath.throughput[2]);
				}
				
			}

			PointKDTree<IptPathState> partialSubPaths(partialSubPathList);

#pragma omp parallel for
			for(int p=0; p<cameraPathNum; p++)
			{
				Path eyePath;
				samplePath(eyePath, camera.generateRay(p));

				if (eyePath.size() <= 1)
					continue;

				IptPathState cameraState;
				cameraState.isSpecularPath = 1;

				cameraState.throughput = vec3f(1.f) / (eyePath[0].directionProb * eyePath[1].originProb);
	
				cameraState.index = eyePath.front().pixelID;

				//cameraState.accuProb = initialProb / (eyePath[0].directionProb);
				cameraState.accuProb = initialProb;

				for(unsigned i=1; i<eyePath.size(); i++)
				{
					vec3f colorHitLight(0.f) , colorDirIllu(0.f) , colorGlbIllu(0.f);

					Real dist = std::max((eyePath[i].origin - eyePath[i - 1].origin).length() , 1e-10f);
					cameraState.throughput *= eyePath[i - 1].getRadianceDecay(dist);

					/*
					cameraState.accuProb *= eyePath[i - 1].directionProb * eyePath[i].originProb;
					cameraState.accuProb /= dist * dist;
					vec3f dir = eyePath[i - 1].origin - eyePath[i].origin;
					dir.normalize();
					cameraState.accuProb *= std::max(0.f , eyePath[i].getContactNormal().dot(dir));
					*/
					/*
					if (cameraState.accuProb > 1e50)
					{
						fprintf(fp , "=========length = %d===========\n", i);
						fprintf(fp , "accuProb = %.10lf , prob[i-1]=%.8f, prob[i]=%.8f, dist=%.8f, cos=%.8f\n" , 
							cameraState.accuProb , eyePath[i - 1].directionProb , eyePath[i].originProb , dist , eyePath[i].getContactNormal().dot(dir));
					}
					*/

					if(eyePath[i].contactObject && eyePath[i].contactObject->emissive())
					{
						if (cameraState.isSpecularPath)
						{
							vec3f le = ((SceneEmissiveObject*)(eyePath[i].contactObject))->getColor();
							colorHitLight = le * cameraState.throughput;
							omp_set_lock(&pixelLocks[cameraState.index]);
							singleImageColors[cameraState.index] += colorHitLight;
							omp_unset_lock(&pixelLocks[cameraState.index]);
						}	
						break;
					}

					cameraState.pos = eyePath[i].origin;
					cameraState.lastRay = &eyePath[i - 1];
					cameraState.ray = &eyePath[i];

					if(eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						colorDirIllu = colorByConnectingLights(camera , singleImageColors , cameraState);
						colorGlbIllu = colorByMergingPaths(singleImageColors , cameraState , partialSubPaths);
					}

					omp_set_lock(&pixelLocks[cameraState.index]);
					singleImageColors[cameraState.index] += colorDirIllu;
					singleImageColors[cameraState.index] += colorGlbIllu;
					if (_isnan(y(colorGlbIllu)))
					{
						printf("indir image error\n");
					}
						
					if (_isnan(y(colorDirIllu)))
					{
						printf("dir image error\n");
					}
					omp_unset_lock(&pixelLocks[cameraState.index]);

					if (eyePath[i].contactObject != NULL && eyePath[i].directionSampleType == Ray::RANDOM)
					{
						cameraState.isSpecularPath = 0;
					}

					if (i >= eyePath.size() - 1)
						break;

					Ray inRay = eyePath[i + 1];
					inRay.direction = -eyePath[i].direction;
					Ray outRay = eyePath[i];
					outRay.direction = -eyePath[i - 1].direction;
				
					Real pdf;
					if(eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						pdf = inRay.getDirectionSampleProbDensity(outRay);
					}
					else
					{
						pdf = eyePath[i].directionProb;
					}
				
					vec3f bsdfFactor;
					bsdfFactor = eyePath[i].color;
				
					cameraState.throughput *= (bsdfFactor *
						eyePath[i].getCosineTerm() / 
						(eyePath[i + 1].originProb * eyePath[i].directionProb));
				
					if (_isnan(y(cameraState.throughput)))
					{
						printf("eye thr error!\n");
					}

					if (eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						Real weightFactor;

						Real volMergeScale = 1;
						if (eyePath[i].insideObject && eyePath[i].contactObject == NULL)
						{
							volMergeScale = 4.0 / 3.0 * mergeRadius;
						}

						weightFactor = connectFactor(pdf) /
							(connectFactor(pdf) + mergeFactor(&volMergeScale , NULL , &INV_2_PI));
						cameraState.throughput *= weightFactor;
					}
				}
			}

			if(cmd == "exit")
				return pixelColors;

			eliminateVignetting(singleImageColors);

			for(int i=0; i<pixelColors.size(); i++)
			{
				pixelColors[i] *= s / Real(s + 1);
				pixelColors[i] += singleImageColors[i] / (s + 1) * camera.width * camera.height; 
			}
			for (int i = 0; i < lightPathNum; i++)
			{
				delete lightPathList[i];
			}
			
			for (int i = 0; i < interPathNum; i++)
			{
				delete interPathList[i];
			}
			
			printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, clock()/1000);

			showCurrentResult(pixelColors);
		}
		else
		{
		}
	}

	for(int i=0; i<pixelLocks.size(); i++)
	{
		omp_destroy_lock(&pixelLocks[i]);
	}
	omp_destroy_lock(&cmdLock);

	return pixelColors;
}

void IptTracer::genLightPaths(omp_lock_t& cmdLock , vector<Path*>& lightPathList)
{
#pragma omp parallel for
	for(int p=0; p<lightPathNum; p++)
	{
		Ray lightRay = genEmissiveSurfaceSample();
		lightPathList[p] = new Path;
		Path &lightPath = *lightPathList[p];
		samplePath(lightPath, lightRay);

		partPathMergeIndex[p].clear();

		if (lightPath.size() <= 1)
			continue;

		IptPathState lightState;
		lightState.isSpecularPath = 1;
		lightState.originRay = &lightPath[0];

		Real cosAtLight = lightPath[0].getCosineTerm();

		lightState.throughput = vec3f(1.0);
		lightState.dirContrib = lightPath[0].color * cosAtLight / 
			(lightPath[0].originProb * lightPath[0].directionProb *
			lightPath[1].originProb);
		lightState.indirContrib = vec3f(0.0);

		for(unsigned i=1; i<lightPath.size(); i++)
		{
			Real dist = std::max((lightPath[i].origin - lightPath[i - 1].origin).length() , 1e-10f);
			vec3f decayFactor = lightPath[i - 1].getRadianceDecay(dist);
			lightState.throughput *= decayFactor;
			lightState.dirContrib *= decayFactor;

			if(lightPath[i].contactObject && lightPath[i].contactObject->emissive())
				break;

			lightState.pos = lightPath[i].origin;
			lightState.lastRay = &lightPath[i - 1];
			lightState.ray = &lightPath[i];

			if(lightPath[i].directionSampleType != Ray::DEFINITE &&
				(lightPath[i].insideObject != NULL || lightPath[i].contactObject != NULL) &&
				(lightState.pos != lightState.originRay->origin))
			{
				omp_set_lock(&cmdLock);
				partialSubPathList.push_back(lightState);
				partPathMergeIndex[p].push_back(partialSubPathList.size() - 1);
				omp_unset_lock(&cmdLock);
			}

			lightState.isSpecularPath &= (lightPath[i].directionSampleType == Ray::DEFINITE);
			Real pdf = lightPath[i].directionProb;

			if (i == lightPath.size() - 1)
				break;

			vec3f scatterFactor = (lightPath[i].color * 
				lightPath[i].getCosineTerm() / 
				(lightPath[i + 1].originProb * lightPath[i].directionProb));

			lightState.throughput *= scatterFactor;
			lightState.dirContrib *= scatterFactor;

			if (_isnan(y(lightState.throughput)))
			{
				printf("light thr error!\n");
			}

			if (lightPath[i].directionSampleType != Ray::DEFINITE)
			{
				Real weightFactor;

				Real volMergeScale = 1;
				if (lightPath[i].insideObject && lightPath[i].contactObject == NULL)
				{
					volMergeScale = 4.0 / 3.0 * mergeRadius;
				}

				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale , NULL , &INV_2_PI));
				lightState.throughput *= weightFactor;
				lightState.dirContrib *= weightFactor;
			}
		}
	}
}

void IptTracer::genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>& interPathList)
{
	// preprocess
	vector<IptPathState> lightSubPathList(partialSubPathList);
	int N = lightSubPathList.size();
	weights.resize(N + 1 , 0);
	for (int i = 0; i < N; i++)
	{
		float intensity = y(lightSubPathList[i].dirContrib);
		weights[i + 1] = weights[i] + intensity * intensity;
	}
	float sum = weights[N];
	for (int i = 0; i <= N; i++)
		weights[i] /= sum;

#pragma omp parallel for
	for(int p=0; p<interPathNum; p++)
	{
		Ray interRay = genIntermediateSamples(lightSubPathList ,
			renderer->scene);

		interPathList[p] = new Path;
		Path &interPath = *interPathList[p];
		samplePath(interPath, interRay);

		partPathMergeIndex[lightPathNum + p].clear();

		if (interPath.size() <= 1)
			continue;

		IptPathState interState;
		interState.isSpecularPath = 1;
		interState.originRay = &interPath[0];

		interState.throughput = vec3f(1.f) / (interPath[0].originProb *
			interPath[0].directionProb * interPath[1].originProb);
		interState.dirContrib = interState.indirContrib = vec3f(0.0);

		for(unsigned i=1; i<interPath.size(); i++)
		{
			Real dist = std::max((interPath[i].origin - interPath[i - 1].origin).length() , 1e-10f);
			interState.throughput *= interPath[i - 1].getRadianceDecay(dist);

			if(interPath[i].contactObject && interPath[i].contactObject->emissive())
				break;

			interState.pos = interPath[i].origin;
			interState.lastRay = &interPath[i - 1];
			interState.ray = &interPath[i];
			/*
			if (y(interState.throughput) > 10)
			{
				fprintf(fp , "throughput = (%.8f,%.8f,%.8f)\n" , interState.throughput[0] , 
					interState.throughput[1] , interState.throughput[2]);
			}
			*/
			if(interPath[i].directionSampleType != Ray::DEFINITE  &&
				(interPath[i].insideObject != NULL || interPath[i].contactObject != NULL) &&
				(interState.pos != interState.originRay->origin))
			{
				if (interPath[i].contactObject && 
					interPath[i].contactObject->emissive())
					break;
				omp_set_lock(&cmdLock);
				partialSubPathList.push_back(interState);
				partPathMergeIndex[lightPathNum + p].push_back(partialSubPathList.size() - 1);
				omp_unset_lock(&cmdLock);
			}

			if (i == interPath.size() - 1)
				break;

			interState.isSpecularPath &= (interPath[i].directionSampleType == Ray::DEFINITE);

			Real pdf = interPath[i].directionProb;

			interState.throughput *= (interPath[i].color * 
				interPath[i].getCosineTerm() / 
				(interPath[i + 1].originProb * interPath[i].directionProb));

			if (_isnan(y(interState.throughput)))
			{
				printf("eye thr error!\n");
			}

			if (interPath[i].directionSampleType != Ray::DEFINITE)
			{
				Real weightFactor;

				Real volMergeScale = 1;
				if (interPath[i].insideObject && interPath[i].contactObject == NULL)
				{
					volMergeScale = 4.0 / 3.0 * mergeRadius;
				}

				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale , NULL , &INV_2_PI));
				interState.throughput *= weightFactor;
			}
		}
	}
}

void IptTracer::mergePartialPaths(omp_lock_t& cmdLock)
{
	struct SearchQuery
	{
		const IptPathState* interState;
		vector<int> mergeIndex;

		SearchQuery() 
		{ 
			mergeIndex.clear();
		}

		void process(const IptPathState& lightState)
		{
			mergeIndex.push_back(lightState.index);
		}
	};

	vector<vec3f> contribs(partialSubPathList.size());

	for (int i = 0; i < partialSubPathList.size(); i++)
	{
		partialSubPathList[i].index = i;
	}

	// preprocess
	PointKDTree<IptPathState> lightTree(partialSubPathList);

	revIndex = new int[partialSubPathList.size()];

	int flag = 0;

#pragma omp parallel for
	for (int i = 0; i < lightPathNum + interPathNum; i++)
	{
		SearchQuery query;

		if (partPathMergeIndex[i].size() == 0)
			continue;

		for (int j = 0; j < partPathMergeIndex[i].size(); j++)
		{
			int k = partPathMergeIndex[i][j];
			revIndex[k] = i;
		}

		query.interState = &partialSubPathList[partPathMergeIndex[i][0]];

		//fprintf(fp , "===========\n");
		//fprintf(fp , "interPos = (%.8f,%.8f,%.8f)\n" , query.interState->pos[0] , 
		//	query.interState->pos[1] , query.interState->pos[2]);

		lightTree.searchInRadius(0 , query.interState->originRay->origin , mergeRadius , query);
		//fprintf(fp , "mergeNum = %d\n" , query.mergeIndex.size());

		partPathMergeIndex[i].clear();
		/*
		omp_set_lock(&cmdLock);
		if (query.mergeIndex.size() > 100)
		{
			flag++;
		}
		*/
		for (int j = 0; j < query.mergeIndex.size(); j++)
		{
			partPathMergeIndex[i].push_back(query.mergeIndex[j]);
			/*
			if (flag == 1 && query.mergeIndex.size() > 100)
			{
			IptPathState& subPath = partialSubPathList[query.mergeIndex[j]];
			
			if (subPath.ray->insideObject && !subPath.ray->contactObject)
			{
				fprintf(fp , "volume\n");
			}
			else if (subPath.ray->contactObject)
			{
				if (subPath.ray->contactObject->emissive())
					fprintf(fp , "light\n");
				else
					fprintf(fp , "surface\n");
			}
			fprintf(fp , "pos=(%.4f,%.4f,%.4f), originPos=(%.4f,%.4f,%.4f)\ndirContrib=(%.8f,%.8f,%.8f)\n" ,
				subPath.pos[0] , subPath.pos[1] , subPath.pos[2] ,
				subPath.originRay->origin[0] , subPath.originRay->origin[1] , subPath.originRay->origin[2] ,
				subPath.dirContrib[0] , subPath.dirContrib[1] , subPath.dirContrib[2]);
			}
			*/
		}
		//omp_unset_lock(&cmdLock);
	}

	int mergeIterations = 10;

	for (int mergeIter = 0; mergeIter < mergeIterations; mergeIter++)
	{
#pragma omp parallel for
		for (int i = 0; i < partialSubPathList.size(); i++)
		{
			mergePartialPaths(contribs , partialSubPathList[i]);
		}

		for (int i = 0; i < partialSubPathList.size(); i++)
		{
			partialSubPathList[i].indirContrib = contribs[i];
		}
	}

	for (int i = 0; i < partPathMergeIndex.size(); i++)
	{
		partPathMergeIndex[i].clear();
		partPathMergeIndex[i].shrink_to_fit();
	}
	partPathMergeIndex.clear();

	delete[] revIndex;
}

Ray IptTracer::genIntermediateSamples(vector<IptPathState>& partialSubPathList , Scene& scene)
{
	float randWeight = RandGenerator::genFloat();
	int pathId = (lower_bound(weights.begin() , weights.end() , randWeight) - weights.begin()) - 1;
	IptPathState& lightState = partialSubPathList[pathId];
	while (lightState.ray == NULL ||
		(lightState.ray->insideObject && !lightState.ray->insideObject->canMerge) ||
		(!lightState.ray->insideObject && lightState.ray->contactObject && !lightState.ray->contactObject->canMerge))
	{
		randWeight = RandGenerator::genFloat();
		pathId = (lower_bound(weights.begin() , weights.end() , randWeight) - weights.begin()) - 1;
		lightState = partialSubPathList[pathId];
	}
	
	Ray ray;
	ray.originSampleType = Ray::SampleType::RANDOM;
	ray.directionSampleType = Ray::SampleType::RANDOM;

	ray.insideObject = lightState.ray->insideObject;
	ray.contactObject = lightState.ray->contactObject;
	ray.contactObjectTriangleID = lightState.ray->contactObjectTriangleID;

	RandGenerator rng;
	Ray inRay , outRay(ray);
	inRay = *lightState.lastRay;

	vec3f o = lightState.pos;
	vec3f dir(0.f);

	if (inRay.insideObject)
	{
		outRay = inRay.insideObject->scatter(inRay);
		dir = outRay.direction;
	}
	else if (inRay.intersectObject)
	{
		outRay = inRay.intersectObject->scatter(inRay);
		dir = outRay.direction;
	}
	
	if (dir.length() < 1e-7f)
	{
		dir = rng.genSphericalDirection();
		outRay.directionProb = 0.25f * M_PI;
	}
	dir.normalize();
	ray.origin = o + (dir * 1e-4 * rng.genFloat());

	ray.direction = dir;

	ray.color = vec3f(1.0);

	ray.originProb = (weights[pathId + 1] - weights[pathId]);

	if (pathId < 0 || pathId >= weights.size())
	{
		printf("pathId error\n");
	}

	//ray.originProb = 1.f;
	ray.directionProb = outRay.directionProb;
	
	ray.current_tid = scene.getContactTreeTid(ray);
	Scene::ObjSourceInformation osi;
	NoSelfIntersectionCondition condition(&scene , ray);
	Real dist = scene.intersect(ray, osi, &condition);

	if (dist > 0)
	{
		ray.intersectDist = dist;
		ray.intersectObject = scene.objects[osi.objID];
		ray.intersectObjectTriangleID = osi.triangleID;
		if (ray.intersectObject == NULL && osi.objID >= 0)
			printf("error\n");
	}

	return ray;
}

void IptTracer::colorByConnectingCamera(vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const IptPathState& lightState)
{
	vec2<float> pCoord = camera.transToPixel(lightState.pos);
	int x = pCoord.x;
	int y = pCoord.y;

	if(!(x >= 0 && x < camera.width && y >= 0 && y < camera.height))
		return;

	vec3f dirToCamera = camera.position - lightState.pos;
	vec3f forward = camera.focus - camera.position;
	forward.normalize();

	if (forward.dot(-dirToCamera) <= 0)
		return;

	Real dist = dirToCamera.length();
	Real distEye2 = dist * dist;
	Real cameraDistToScreen2 = camera.sightDist * camera.sightDist;
	dirToCamera = dirToCamera / dist;

	Real cosToCamera;
	if (lightState.ray->insideObject && !lightState.ray->contactObject)
	{
		cosToCamera = 1.f;
	}
	else
	{
		cosToCamera = std::abs(lightState.ray->getContactNormal().dot(dirToCamera));
	}

	Ray outRay = *lightState.ray;
	outRay.direction = dirToCamera;

	vec3f decayFactor = outRay.getRadianceDecay(dist);

	vec3f bsdfFactor = lightState.lastRay->getBSDF(outRay);

	Real cosAtCamera = forward.dot(-dirToCamera);

	Real volMergeScale = 1;
	if (lightState.ray->insideObject && lightState.ray->contactObject == NULL)
	{
		volMergeScale = 4.0 / 3.0 * mergeRadius;
	}

	Real imagePointToCameraDist = camera.sightDist / cosAtCamera;
	Real imageToSolidAngleFactor = imagePointToCameraDist *
		imagePointToCameraDist / cosAtCamera;
	Real imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / distEye2;

	Real cameraPdfArea = imageToSurfaceFactor * 1.f; // pixel area is 1
	
	Real surfaceToImageFactor = 1.f / imageToSurfaceFactor;
	
	Real bsdfDirPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);

	vec3f totContrib = lightState.dirContrib + lightState.indirContrib;

	//---- still buggy, fix me ----
	vec3f color = (totContrib * bsdfFactor * decayFactor) /
		(lightPathNum * surfaceToImageFactor);
	color *= powf(cosAtCamera , 4) / lightPathNum;

	//vec3f color = (totContrib * bsdfFactor * decayFactor) * 
	//	cosAtCamera * cosToCamera * cameraDistToScreen2 / 
	//	(distEye2 * lightPathNum * lightPathNum);

	//vec3f color = (totContrib * bsdfFactor * decayFactor) * 
	//	cosAtCamera * cosToCamera / (cameraDistToScreen2 * lightPathNum);
	//-----------------------------

	Ray inRay;
	inRay.origin = camera.position;
	inRay.direction = -dirToCamera;

	if (!testVisibility(inRay , outRay))
		return;

	Real pdf = bsdfDirPdf;
	
	Real weightFactor = connectFactor(pdf) / 
		(connectFactor(pdf) + mergeFactor(&volMergeScale));

	//color *= weightFactor;

// 	fprintf(fp , "weight = %.6lf, connect camera = (%.6lf,%.6lf,%.6lf)\n" ,
// 		weightFactor , color[0] , color[1] , color[2]);

	omp_set_lock(&pixelLocks[y*camera.width + x]);
	colors[y*camera.width + x] += color * 0.5f;
	omp_unset_lock(&pixelLocks[y*camera.width + x]);
}

vec3f IptTracer::colorByConnectingLights(const Camera& camera, vector<vec3f>& colors, const IptPathState& cameraState)
{
	Ray lightRay = genEmissiveSurfaceSample();
	lightRay.direction = (cameraState.pos - lightRay.origin);
	Real dist2 = lightRay.direction.length();
	dist2 = dist2 * dist2;
	lightRay.direction.normalize();
	Ray outRay = *cameraState.ray;
	outRay.direction = -lightRay.direction;

	vec3f decayFactor = outRay.getRadianceDecay(sqrt(std::max(0.f , dist2)));

	if(!testVisibility(outRay, lightRay))
		return vec3f(0.f);

	//outRay.direction = -cameraState.lastRay->direction;
	//vec3f bsdfFactor = lightRay.getBSDF(outRay);
	vec3f bsdfFactor = cameraState.lastRay->getBSDF(outRay);
	
	Real cosAtLight = min(max(0.f , lightRay.getContactNormal().dot(lightRay.direction)) , 1.f);

	Real cosToLight = 1;
	if (cameraState.ray->contactObject != NULL)
	{
		cosToLight = min(max(0.f , cameraState.ray->getContactNormal().dot(-lightRay.direction)) , 1.f);
	}

	Real volMergeScale = 1;
	if (cameraState.ray->insideObject && cameraState.ray->contactObject == NULL)
	{
		volMergeScale = 4.0 / 3.0 * mergeRadius;
	}

	vec3f tmp = lightRay.color * cosAtLight * bsdfFactor * cosToLight * cameraState.throughput
		/ (lightRay.originProb * dist2);

	outRay.direction = -cameraState.lastRay->direction;
	Real pdf = lightRay.getDirectionSampleProbDensity(outRay);

	Real weightFactor = connectFactor(pdf) /
		(connectFactor(pdf) + mergeFactor(&volMergeScale , NULL , &INV_2_PI));

	vec3f res = tmp * decayFactor * weightFactor;
	/*
	vec3f resx = camera.eliminateVignetting(res , cameraState.index) * lightPathNum;
	if (resx[0] + resx[1] + resx[2] > 0)
	{
		fprintf(fp , "=====================\n");
		fprintf(fp , "decay=(%.4f,%.4f,%.4f)\nscatter=(%.4f,%.4f,%.4f)\nlight=(%.4f,%.4f,%.4f)\nweight=%.4f, res=(%.10f,%.10f,%.10f)\n" , 
			decayFactor[0] , decayFactor[1] , decayFactor[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] ,
			lightRay.color[0] , lightRay.color[1] , lightRay.color[2] , weightFactor , resx[0] , resx[1] , resx[2]);
	}
	*/
	return res;
}

vec3f IptTracer::colorByMergingPaths(vector<vec3f>& colors, const IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths)
{
	struct GatherQuery
	{
		vec3f color;
		IptTracer *tracer;
		const IptPathState* cameraState;
		int mergeNum;

		GatherQuery(IptTracer* tracer) { this->tracer = tracer; mergeNum = 0; }

		void process(const IptPathState& lightState)
		{
			Real volMergeScale = 1;

			if (lightState.ray->insideObject && lightState.ray->contactObject == NULL)
			{
				if (cameraState->ray->insideObject == NULL || 
					cameraState->ray->contactObject != NULL ||
					cameraState->ray->insideObject != lightState.ray->insideObject ||
					!cameraState->ray->insideObject->canMerge ||
					!lightState.ray->insideObject->canMerge)
				{
					return;
				}
				volMergeScale = 4.0 / 3.0 * tracer->mergeRadius;
			}
			else if (lightState.ray->contactObject)
			{
				if (lightState.ray->contactObject != cameraState->ray->contactObject ||
					!lightState.ray->contactObject->canMerge ||
					!cameraState->ray->contactObject->canMerge)
				{
					return;
				}
			}
			else 
			{
				return;
			}

			Ray outRay;
			vec3f bsdfFactor;
			
			outRay = *lightState.ray;
			outRay.direction = -(cameraState->lastRay->direction);
			bsdfFactor = lightState.lastRay->getBSDF(outRay);
			if (y(bsdfFactor) < 1e-8)
				return;

			vec3f totContrib = lightState.dirContrib + lightState.indirContrib; 
			vec3f tmp = totContrib * bsdfFactor * 
				cameraState->throughput; 

			Real lastPdf , weightFactor;
			lastPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);

			Real originProb = min(cameraState->accuProb , 1e10);
			weightFactor = tracer->mergeFactor(&volMergeScale , &originProb , &INV_2_PI) / 
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &originProb , &INV_2_PI));

			mergeNum++;
			vec3f res = tmp * (weightFactor * tracer->mergeKernel / volMergeScale);
			color += res;
			
			vec3f resx = tracer->renderer->camera.eliminateVignetting(res , cameraState->index) *
				tracer->pixelNum;	
			
			if (_isnan(y(res)))
			{
				fprintf(fp , "===========nan==========\n");
				if (volMergeScale == 1)
					fprintf(fp , "surface\n");
				else 
					fprintf(fp , "volume\n");
				fprintf(fp , "res = (%.8f,%.8f,%.8f) \ntotContrib = (%.8f,%.8f,%.8f), bsdf = (%.8f,%.8f,%.8f),\n cameraThr = (%.8f,%.8f,%.8f) \nweightFactor = %.8f, originProb = %.8f, lastpdf = %.8f\n" ,
					resx[0] , resx[1] , resx[2] ,
					totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , 
					cameraState->throughput[0] , cameraState->throughput[1] , cameraState->throughput[2] , 
					weightFactor , cameraState->accuProb , lastPdf);
			}
			
		}
	};

	GatherQuery query(this);
	query.cameraState = &cameraState;
	query.color = vec3f(0, 0, 0);

	partialSubPaths.searchInRadius(0 , query.cameraState->pos , mergeRadius , query);
	
	//fprintf(fp , "mergeNum = %d\n" , query.mergeNum);

	return query.color;
}


void IptTracer::mergePartialPaths(vector<vec3f>& contribs , const IptPathState& lightState)
{
	struct MergeQuery
	{
		vec3f color;
		IptTracer *tracer;
		const IptPathState* interState;

		MergeQuery(IptTracer* tracer) { this->tracer = tracer; }

		void process(const IptPathState& lightState)
		{
			Real volMergeScale = 1;

			if (lightState.ray->insideObject && lightState.ray->contactObject == NULL)
			{
				if (interState->originRay->insideObject == NULL || 
					interState->originRay->contactObject != NULL ||
					interState->originRay->insideObject != lightState.ray->insideObject ||
					!interState->originRay->insideObject->canMerge ||
					!lightState.ray->insideObject->canMerge)
				{
					return;
				}
				volMergeScale = 4.0 / 3.0 * tracer->mergeRadius;
			}
			else if (lightState.ray->contactObject)
			{
				if (lightState.ray->contactObject != interState->originRay->contactObject ||
					!lightState.ray->contactObject->canMerge ||
					!interState->originRay->contactObject->canMerge)
				{
					return;
				}
			}
			else
			{
				return;
			}
			
			vec3f totContrib = lightState.dirContrib + lightState.indirContrib;

			Ray outRay;
			vec3f bsdfFactor;
		
			outRay = *lightState.ray;
			outRay.direction = interState->originRay->direction;
			bsdfFactor = lightState.lastRay->getBSDF(outRay);

			vec3f tmp = totContrib * bsdfFactor * 
				interState->throughput;

			Real lastPdf , weightFactor;
			lastPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);
			weightFactor = tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb ,
				&interState->originRay->directionProb) / 
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb ,
				&interState->originRay->directionProb));
			
			color += tmp * (weightFactor * tracer->mergeKernel / volMergeScale);
			/*
			vec3f resx = tmp * weightFactor;	
			if (resx[0] + resx[1] + resx[2] > 10)
			{
				fprintf(fp , "=====================\n");
				fprintf(fp , "res = (%.8f,%.8f,%.8f), totContrib = (%.8f,%.8f,%.8f), \nbsdf = (%.8f,%.8f,%.8f), \ninterThr = (%.8f,%.8f,%.8f), weightFactor = %.8f\n" ,
					resx[0] , resx[1] , resx[2] , totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , 
					interState->throughput[0] , interState->throughput[1] , interState->throughput[2] , weightFactor);
			}
			*/
		}
	};

	MergeQuery query(this);
	query.interState = &lightState;
	query.color = vec3f(0, 0, 0);

	int pa = revIndex[lightState.index];
	//fprintf(fp , "%d\n" , partPathMergeIndex[pa].size());
	for (int j = 0; j < partPathMergeIndex[pa].size(); j++)
	{
		int k = partPathMergeIndex[pa][j];
		query.process(partialSubPathList[k]);
	}

	contribs[lightState.index] = query.color;
}