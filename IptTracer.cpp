#include "StdAfx.h"
#include "IptTracer.h"
#include "SceneEmissiveObject.h"
#include "UniformSphericalSampler.h"
#include "NoSelfIntersectionCondition.h"

static FILE *fp = fopen("debug_ipt.txt" , "w");

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

			mergeRadius = r0 * sqrt(pow(s+1, alpha-1));
			mergeRadius = std::max(mergeRadius , 1e-7f);

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

			printf("%d\n" , partialSubPathList.size());

#pragma omp parallel for
			for (int i = 0; i < partialSubPathList.size(); i++)
			{
				IptPathState& subPath = partialSubPathList[i];
				colorByConnectingCamera(pixelLocks , camera , singleImageColors , subPath);

				/*
				if (s == 0)
				{
				fprintf(fp , "dirContrib=(%.4f,%.4f,%.4f), indirContrib=(%.4f,%.4f,%.4f), throughput=(%.4f,%.4f,%.4f)\n" ,
					subPath.dirContrib[0] , subPath.dirContrib[1] , subPath.dirContrib[2] ,
					subPath.indirContrib[0] , subPath.indirContrib[1] , subPath.indirContrib[2] ,
					subPath.throughput[0] , subPath.throughput[1] , subPath.throughput[2]);
				}
				*/
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

				//cameraState.throughput = vec3f(1.f) / (eyePath[0].directionProb * eyePath[1].originProb);
				cameraState.throughput = vec3f(1.f) * powf(eyePath[0].getCosineTerm() , 4)
					/ (pixelNum * eyePath[1].originProb);
	
				cameraState.index = eyePath.front().pixelID;

				for(unsigned i=1; i<eyePath.size(); i++)
				{
					vec3f colorHitLight(0.f) , colorDirIllu(0.f) , colorGlbIllu(0.f);

					Real dist = std::max((eyePath[i].origin - eyePath[i - 1].origin).length() , 1e-10f);
					cameraState.throughput *= eyePath[i - 1].getRadianceDecay(dist);

					if(eyePath[i].contactObject && eyePath[i].contactObject->emissive())
					{
						if (cameraState.isSpecularPath)
						{
							vec3f le = ((SceneEmissiveObject*)(eyePath[i].contactObject))->getColor();
							//vec3f dirAtLight = eyePath[i - 1].origin - eyePath[i].origin;
							//dirAtLight.normalize();
							//Real cosAtLight = eyePath[i].getContactNormal().dot(dirAtLight);
							//cosAtLight = clampf(cosAtLight , 0.f , 1.f);

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
					singleImageColors[cameraState.index] += colorDirIllu + colorGlbIllu;
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
				
					if (eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						Real weightFactor;

						Real volMergeScale = 1;
						if (eyePath[i].insideObject && eyePath[i].contactObject == NULL)
						{
							volMergeScale = 4.0 / 3.0 * mergeRadius;
						}

						weightFactor = connectFactor(pdf) /
							(connectFactor(pdf) + mergeFactor(&volMergeScale));
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
				lightPathList[i]->clear();
				delete lightPathList[i];
			}
			
			for (int i = 0; i < interPathNum; i++)
			{
				interPathList[i]->clear();
				delete interPathList[i];
			}
			
			lightPathList.clear();
			interPathList.clear();
			
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

			if (lightPath[i].directionSampleType != Ray::DEFINITE)
			{
				Real weightFactor;

				Real volMergeScale = 1;
				if (lightPath[i].insideObject && lightPath[i].contactObject == NULL)
				{
					volMergeScale = 4.0 / 3.0 * mergeRadius;
				}

				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale));
				lightState.throughput *= weightFactor;
				lightState.dirContrib *= weightFactor;
			}
		}
	}
}

void IptTracer::genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>& interPathList)
{
	vector<IptPathState> lightSubPathList(partialSubPathList);
#pragma omp parallel for
	for(int p=0; p<interPathNum; p++)
	{
		Ray interRay = genIntermediateSamples(lightSubPathList ,
			renderer->scene);

		//Ray interRay = genOtherSurfaceSample();

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

			if(interPath[i].directionSampleType != Ray::DEFINITE  &&
				(interPath[i].insideObject != NULL || interPath[i].contactObject != NULL) &&
				(interState.pos != interState.originRay->origin))
			{
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

			if (interPath[i].directionSampleType != Ray::DEFINITE)
			{
				Real weightFactor;

				Real volMergeScale = 1;
				if (interPath[i].insideObject && interPath[i].contactObject == NULL)
				{
					volMergeScale = 4.0 / 3.0 * mergeRadius;
				}

				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale));
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

#pragma omp parallel for
	for (int i = 0; i < lightPathNum + interPathNum; i++)
	{
		partPathMergeIndex[i].clear();
		
		SearchQuery query;

		if (partPathMergeIndex[i].size() == 0)
			continue;

		for (int j = 0; j < partPathMergeIndex[i].size(); j++)
		{
			int k = partPathMergeIndex[i][j];
			revIndex[k] = i;
		}

		query.interState = &partialSubPathList[partPathMergeIndex[i][0]];

		lightTree.searchInRadius(0 , query.interState->originRay->origin , mergeRadius , query);

		partPathMergeIndex[i].clear();
		for (int j = 0; j < query.mergeIndex.size(); j++)
		{
			partPathMergeIndex[i].push_back(query.mergeIndex[j]);
			/*
			IptPathState& subPath = partialSubPathList[query.mergeIndex[j]];
			fprintf(fp , "pos=(%.4f,%.4f,%.4f), originPos=(%.4f,%.4f,%.4f)\n" ,
				subPath.pos[0] , subPath.pos[1] , subPath.pos[2] ,
				subPath.originRay->origin[0] , subPath.originRay->origin[1] , subPath.originRay->origin[2]);
			*/
		}
	}

	int mergeIterations = 1;

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

	delete[] revIndex;
}

Ray IptTracer::genIntermediateSamples(vector<IptPathState>& partialSubPathList , Scene& scene)
{
	int N = partialSubPathList.size();
	
	int pathId = (int)RandGenerator::genFloat() * N;
	IptPathState& lightState = partialSubPathList[pathId];
	while (lightState.ray == NULL)
	{
		pathId = (int)RandGenerator::genFloat() * N;
		lightState = partialSubPathList[pathId];
	}
	
	Ray ray;
	ray.originSampleType = Ray::SampleType::RANDOM;
	ray.directionSampleType = Ray::SampleType::RANDOM;

	ray.insideObject = lightState.ray->insideObject;
	ray.contactObject = lightState.ray->contactObject;
	ray.contactObjectTriangleID = lightState.ray->contactObjectTriangleID;

	RandGenerator rng;

	vec3f o = lightState.pos;
	vec3f n;
	if (lightState.ray->contactObject != NULL)
	{
		n = lightState.ray->getContactNormal();
	}
	else
	{
		n = rng.genSphericalDirection();
	}
	n.normalize();
	ray.origin = o + (n * 1e-4);

	vec3f dir = rng.genSphericalDirection();
	if (dir.dot(n) <= 0)
		dir = -dir;
	dir.normalize();
	ray.direction = dir;

	ray.color = vec3f(1.0);
	ray.originProb = 1.0;
	ray.directionProb = 1.f / (2.f * M_PI);

	Scene::ObjSourceInformation osi;
	NoSelfIntersectionCondition condition(&scene , ray);
	Real dist = scene.intersect(ray, osi, &condition);

	if (dist > 0)
	{
		ray.intersectDist = dist;
		ray.intersectObject = scene.objects[osi.objID];
		ray.intersectObjectTriangleID = osi.triangleID;
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
	//	cosAtCamera * cosToCamera / (distEye2 * lightPathNum);
	//-----------------------------

	Ray inRay;
	inRay.origin = camera.position;
	inRay.direction = -dirToCamera;

	if (!testVisibility(inRay , outRay))
		return;

	Real pdf = bsdfDirPdf;
	
	Real weightFactor = connectFactor(pdf) / 
		(connectFactor(pdf) + mergeFactor(&volMergeScale));

	color *= weightFactor;

// 	fprintf(fp , "weight = %.6lf, connect camera = (%.6lf,%.6lf,%.6lf)\n" ,
// 		weightFactor , color[0] , color[1] , color[2]);

	omp_set_lock(&pixelLocks[y*camera.width + x]);
	colors[y*camera.width + x] += color;
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

	vec3f decayFactor = outRay.getRadianceDecay(sqrt(dist2));

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
		(connectFactor(pdf) + mergeFactor(&volMergeScale));

	vec3f res = tmp * decayFactor * weightFactor;
	vec3f resx = camera.eliminateVignetting(res , cameraState.index) * lightPathNum;
	/*
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
			SceneObject *obj = NULL;

			if (lightState.ray->insideObject && lightState.ray->contactObject == NULL)
			{
				if (cameraState->ray->insideObject == NULL || 
					cameraState->ray->contactObject != NULL ||
					cameraState->ray->insideObject != lightState.ray->insideObject)
				{
					return;
				}
				volMergeScale = 4.0 / 3.0 * tracer->mergeRadius;
				obj = lightState.ray->insideObject;
			}
			else
			{
				if (lightState.ray->contactObject != cameraState->ray->contactObject)
				{
					return;
				}
				obj = lightState.ray->contactObject;
			}

			Ray outRay;
			vec3f bsdfFactor;
			
			outRay = *lightState.ray;
			outRay.direction = -(cameraState->lastRay->direction);
			bsdfFactor = lightState.lastRay->getBSDF(outRay);

			Real contProb = 1.f;
			if (obj)
			{
				contProb = obj->getContinueProbability(*lightState.lastRay , outRay);
			}

			vec3f totContrib = lightState.dirContrib + lightState.indirContrib;
			vec3f tmp = totContrib * bsdfFactor * 
				cameraState->throughput * cameraState->ray->originProb / volMergeScale; 

			Real lastPdf , weightFactor;
			lastPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);
			weightFactor = tracer->mergeFactor(&volMergeScale) / 
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale));
			weightFactor *= tracer->mergeKernel / volMergeScale;

			mergeNum++;
			color += tmp * weightFactor;
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
					interState->originRay->insideObject != lightState.ray->insideObject)
				{
					return;
				}
				volMergeScale = 4.0 / 3.0 * tracer->mergeRadius;
			}
			else
			{
				if (lightState.ray->contactObject != interState->originRay->contactObject)
				{
					return;
				}
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
			weightFactor = tracer->mergeFactor(&volMergeScale) / 
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale));
			weightFactor *= tracer->mergeKernel / volMergeScale;

			color += tmp * weightFactor;
		}
	};

	MergeQuery query(this);
	query.interState = &lightState;
	query.color = vec3f(0, 0, 0);

	int pa = revIndex[lightState.index];
	for (int j = 0; j < partPathMergeIndex[pa].size(); j++)
	{
		int k = partPathMergeIndex[pa][j];
		query.process(partialSubPathList[k]);
	}

	contribs[lightState.index] = query.color;
}