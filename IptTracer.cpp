#include "StdAfx.h"
#include "IptTracer.h"
#include "SceneEmissiveObject.h"
#include "UniformSphericalSampler.h"
#include "NoSelfIntersectionCondition.h"
#include "SceneVPMObject.h"
#include "SceneHeteroGeneousVolume.h"

static FILE *fp = fopen("debug_ipt_y.txt" , "w");

float INV_2_PI = 0.5 / M_PI;

vector<vec3f> IptTracer::renderPixels(const Camera& camera)
{
	vector<vec3f> pixelColors(camera.width * camera.height, vec3f(0, 0, 0));
	vector<omp_lock_t> pixelLocks(pixelColors.size());

	preprocessEmissionSampler();
	preprocessOtherSampler();
	preprocessVolumeSampler();
	
	for(int i=0; i<pixelLocks.size(); i++)
	{
		omp_init_lock(&pixelLocks[i]);
	}

	omp_lock_t cmdLock;
	omp_init_lock(&cmdLock);

	Real r0 = mergeRadius;

	totArea = renderer->scene.getTotalArea();
	totVol = renderer->scene.getTotalVolume();
	printf("scene: totArea = %.8f, totVol = %.8f\n" , totArea , totVol);

	for(unsigned s=0; s<spp; s++)
	{
		partPathMergeIndex.resize(interPathNum);

		partialSubPathList.clear();

		float shrinkRatio;
		if (totVol > 1e-6f)
			shrinkRatio = powf(((float)s + alpha) / ((float)s + 1.f) , 1.f / 3.f);
		else
			shrinkRatio = powf(((float)s + alpha) / ((float)s + 1.f) , 1.f / 2.f);
		if (s > 0)
			mergeRadius *= shrinkRatio;

		mergeRadius = std::max(mergeRadius , 1e-7f);
		printf("mergeRadius = %.8f\n" , mergeRadius);

		vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));

		string cmd;

		unsigned t = clock();

		vector<Path*> lightPathList(lightPathNum , NULL);
		vector<Path*> interPathList(interPathNum, NULL);

		if (!renderer->scene.usingGPU())
		{
			genLightPaths(cmdLock , lightPathList);

			if (!usePPM)
				genIntermediatePaths(cmdLock , interPathList);
			
			printf("lightPhotonNum = %d, partialPhotonNum = %d\n" , lightPhotonNum , partialPhotonNum);

			mergeKernel = 1.f / (M_PI * mergeRadius * 
				mergeRadius * (Real)partialPathNum);

			mergePartialPaths(cmdLock);

			PointKDTree<IptPathState> partialSubPaths(partialSubPathList);
			
			/*
			for (int i = 0; i < partialPhotonNum; i++)
			{
				IptPathState& subPath = partialSubPathList[i];
				if (s == 0)
				{
					fprintf(fp , "==============\n");
					vec3f contrib;
					if (i < lightPhotonNum)
					{
						contrib = subPath.throughput;
						fprintf(fp , "dirContrib=(%.8f,%.8f,%.8f)\n" , contrib.x , contrib.y , contrib.z);
					}
					else
					{
						contrib = subPath.indirContrib;
						fprintf(fp , "indirContrib=(%.8f,%.8f,%.8f)\n" , contrib.x , contrib.y , contrib.z);
					}			
				}
			}
			*/

#pragma omp parallel for
			for(int p=0; p<cameraPathNum; p++)
			{
				Path eyePath;
				//samplePath(eyePath, camera.generateRay(p));
				sampleMergePath(eyePath , camera.generateRay(p) , 0);
				singleImageColors[p] += colorByRayMarching(eyePath , partialSubPaths);
				
				// abandon all the rest!
				/*
				if (eyePath.size() <= 1)
					continue;

				IptPathState cameraState;
				cameraState.isSpecularPath = 1;

				cameraState.throughput = vec3f(1.f) / (eyePath[0].originProb * eyePath[0].directionProb * eyePath[1].originProb);
	
				cameraState.index = eyePath.front().pixelID;

				//cameraState.accuProb = 1.f;

				vector<double> probDir , probRev;
				vector<float> weights;
				//calcEyePathProbs(eyePath , probDir , probRev);

				//fprintf(fp , "===================\n");
				//for (int i = 0; i < eyePath.size(); i++)
				//{
				//	fprintf(fp , "l=%d, bsdf=(%.8f,%.8f,%.8f), originPdf=%.8f, dirPdf=%.8f\n" , i , eyePath[i].color.x ,
				//		eyePath[i].color.y , eyePath[i].color.z , eyePath[i].originProb , eyePath[i].directionProb);
				//}

				int nonSpecLength = 0;
				vector<vec3f> mergeContribs;
				vec3f mergeContrib(0.f);

				mergeContribs.clear();
				weights.clear();

				float weightFactor = 1.f;

				for(unsigned i = 1; i < eyePath.size(); i++)
				//for (unsigned i = 1; i < 2; i++)
				{
					vec3f colorHitLight(0.f) , colorDirIllu(0.f) , colorGlbIllu(0.f);

					//Real dist = std::max((eyePath[i].origin - eyePath[i - 1].origin).length() , 1e-5f);
					//cameraState.throughput *= eyePath[i - 1].getRadianceDecay(dist);

					if(eyePath[i].contactObject && eyePath[i].contactObject->emissive())
					{
						if (cameraState.isSpecularPath)
						{
							if (abs(eyePath[i].getContactNormal().dot(-eyePath[i - 1].direction)) > 1e-6f)
							{
								vec3f le = ((SceneEmissiveObject*)(eyePath[i].contactObject))->getColor();
								colorHitLight = le * cameraState.throughput;
								singleImageColors[cameraState.index] += colorHitLight;
							}
						}	
						break;
					}

					if ((eyePath[i].insideObject == NULL && eyePath[i].contactObject == NULL) ||
						(eyePath[i].origin == eyePath[i - 1].origin))
						break;

					cameraState.pos = eyePath[i].origin;
					cameraState.lastRay = &eyePath[i - 1];
					cameraState.ray = &eyePath[i];

					Real eyeWeight;
					//colorDirIllu = colorByConnectingLights(camera , cameraState);
					//singleImageColors[cameraState.index] += colorDirIllu;

					colorGlbIllu = colorByMergingPaths(cameraState , partialSubPaths);
					mergeContribs.push_back(colorGlbIllu);
					weights.push_back(weightFactor);
			
					//weights.push_back(probDir[i] * probRev[i]);
					//fprintf(fp , "%.8f " , weightFactor);

					if (eyePath[i].directionSampleType == Ray::RANDOM)
					{
						cameraState.isSpecularPath = 0;
						nonSpecLength++;
						
						//if (usePPM)
						//if (nonSpecLength == 1)
						//	break; // PPM, eye path length is 1
					}

					if (i >= eyePath.size() - 1)
						break;
				
					vec3f bsdfFactor;
					bsdfFactor = eyePath[i].color;
				
					//cameraState.throughput *= (bsdfFactor * eyePath[i].getCosineTerm() / 
					//	(eyePath[i + 1].originProb * eyePath[i].directionProb));

					//fprintf(fp , "l=%d, thr=(%.8f,%.8f,%.8f), bsdf=(%.8f,%.8f,%.8f), cos=%.8f, prob=%.8f\n" , 
					//	i , cameraState.throughput[0] , cameraState.throughput[1] , cameraState.throughput[2] ,
					//	bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , eyePath[i].getCosineTerm() , eyePath[i].directionProb);
					
					if (eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						Ray inRay = eyePath[i + 1];
						inRay.direction = -eyePath[i].direction;
						Ray outRay = eyePath[i];
						outRay.direction = -eyePath[i - 1].direction;

						Real dirPdf = eyePath[i].directionProb;
						Real revPdf = inRay.getDirectionSampleProbDensity(outRay);
						//printf("%.8f, %.8f\n" , dirPdf , revPdf);

						Real volMergeScale = 1;
						Real originProb = 1.f / totArea;
						Real dirProb = eyePath[i].getCosineTerm() / M_PI;
						if (eyePath[i].insideObject && !eyePath[i].contactObject)
						{
							volMergeScale = 4.0 / 3.0 * mergeRadius;
							originProb = 1.f / totVol;
							dirProb = 0.25f / M_PI;
						}

						weightFactor *= connectFactor(revPdf) /
							(connectFactor(revPdf) + mergeFactor(&volMergeScale , &originProb , &dirProb));

						if (_isnan(weightFactor))
							printf("error, w = %.8f\n" , weightFactor);

						//cameraState.throughput *= weightFactor;
						//weights.push_back(weightFactor);
					}
				}
				
				double sumWeights = 0.f;
				for (int i = 0; i < weights.size(); i++)
					sumWeights += weights[i];
				for (int i = 0; i < mergeContribs.size(); i++)
				{
					singleImageColors[cameraState.index] += mergeContribs[i];// * weights[i] / sumWeights;
					//if (nonSpecLength > 0)
					//	singleImageColors[cameraState.index] += mergeContribs[i] / (Real)nonSpecLength;		
				}
				*/
			}
		}
		else
		{
			vector<Path> lightPathListGPU , interPathListGPU , eyePathListGPU;
			vector<Ray> eyeRayList(cameraPathNum);
			vector<Ray> lightRayList(lightPathNum);
			vector<Ray> interRayList(interPathNum);

#pragma omp parallel for
			for (int p = 0; p < lightPathNum; p++)
				lightRayList[p] = genEmissiveSurfaceSample();
			lightPathListGPU = samplePathList(lightRayList);
			movePaths(cmdLock , lightPathListGPU , lightPathList);

			genLightPaths(cmdLock , lightPathList);

			if (!usePPM)
			{
				vector<IptPathState> lightSubPathList(partialSubPathList);
#pragma omp parallel for
				for (int p = 0; p < interPathNum; p++)
					interRayList[p] = genIntermediateSamples(renderer->scene);
				interPathListGPU = samplePathList(interRayList);
				movePaths(cmdLock , interPathListGPU , interPathList);

				genIntermediatePaths(cmdLock , interPathList);
			}
			
			printf("lightPhotonNum = %d, partialPhotonNum = %d\n" , lightPhotonNum , partialPhotonNum);

			mergeKernel = 1.f / (M_PI * mergeRadius * 
				mergeRadius * (Real)partialPathNum);

			mergePartialPaths(cmdLock);

			PointKDTree<IptPathState> partialSubPaths(partialSubPathList);

//#pragma omp parallel for
//			for (int p = 0; p < cameraPathNum; p++)
//				eyeRayList[p] = camera.generateRay(p);
//			eyePathListGPU = samplePathList(eyeRayList);

#pragma omp parallel for
			for(int p=0; p<cameraPathNum; p++)
			{
				Path eyePath;
				sampleMergePath(eyePath , camera.generateRay(p) , 0);
				singleImageColors[p] += colorByRayMarching(eyePath , partialSubPaths);
			}
		}

		printf("done calculation, release memory\n");

		if(cmd == "exit")
			return pixelColors;

		for(int i=0; i<pixelColors.size(); i++)
		{
			pixelColors[i] *= (Real)s / ((Real)s + 1.f);
			pixelColors[i] += singleImageColors[i] / ((Real)s + 1.f); 
		}

		if (!renderer->scene.usingGPU())
		{
			for (int i = 0; i < lightPathNum; i++)
			{
				if (lightPathList[i])
					delete lightPathList[i];
			}

			for (int i = 0; i < interPathNum; i++)
			{
				if (interPathList[i])
					delete interPathList[i];
			}
		}
		else
		{
			for (int i = 0; i < lightPathNum; i++)
				lightPathList[i] = NULL;
			for (int i = 0; i < interPathNum; i++)
				interPathList[i] = NULL;
		}

		printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, clock()/1000);

		if (clock() / 1000 >= lastTime)
		{
			showCurrentResult(pixelColors , &lastTime);
			lastTime += timeInterval;
		}
		else
			showCurrentResult(pixelColors);
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
		if (!renderer->scene.usingGPU())
		{
			Ray lightRay = genEmissiveSurfaceSample();
			lightPathList[p] = new Path;
			samplePath(*lightPathList[p] , lightRay);
		}
		
		Path& lightPath = *lightPathList[p];

		if (lightPath.size() <= 1)
			continue;

		IptPathState lightState;
		lightState.originRay = &lightPath[0];

		Real cosAtLight = lightPath[0].getCosineTerm();

		lightState.throughput = lightPath[0].color * cosAtLight / 
			(lightPath[0].originProb * lightPath[0].directionProb *
			lightPath[1].originProb);
		lightState.indirContrib = vec3f(0.0);
		/*
		fprintf(fp , "====================\n");
		vec3f decay = lightPath[0].getRadianceDecay((lightPath[0].origin - lightPath[1].origin).length());
		fprintf(fp , "l = 0 , thr = (%.8f,%.8f,%.8f) , color = (%.8f,%.8f,%.8f)\ncosine = %.8f , dirPdf = %.8f , oPdf = %.8f\ndecay=(%.8f,%.8f,%.8f)\n" , 
			lightState.throughput.x , lightState.throughput.y , lightState.throughput.z ,
			lightPath[0].color[0] , lightPath[0].color[1] , lightPath[0].color[2] ,
			lightPath[0].getCosineTerm() , lightPath[0].directionProb , lightPath[1].originProb ,
			decay.x , decay.y , decay.z);
		*/
		int nonSpecPathLength = 0;

		for(unsigned i = 1; i < lightPath.size(); i++)
		//for (unsigned i = 1; i < 2; i++)
		{
			Real dist = std::max((lightPath[i].origin - lightPath[i - 1].origin).length() , 1e-5f);
			vec3f decayFactor = lightPath[i - 1].getRadianceDecay(dist);
			lightState.throughput *= decayFactor;

			if(lightPath[i].contactObject && lightPath[i].contactObject->emissive())
				break;

			lightState.pos = lightPath[i].origin;
			lightState.lastRay = &lightPath[i - 1];
			lightState.ray = &lightPath[i];
			lightState.pathLen = i;

			if(lightPath[i].directionSampleType == Ray::RANDOM &&
				(lightPath[i].insideObject != NULL || lightPath[i].contactObject != NULL) && 
				(lightPath[i].origin != lightPath[i - 1].origin))
			{
				//if (lightPath[i].insideObject && !lightPath[i].contactObject)
				//	fprintf(fp , "path length = %d, dirContrib = (%.8f,%.8f,%.8f)\n" , 
				//		i , lightState.dirContrib[0] , lightState.dirContrib[1] , lightState.dirContrib[2]);

				omp_set_lock(&cmdLock);
				partialSubPathList.push_back(lightState);
				omp_unset_lock(&cmdLock);
			}

			if (i == lightPath.size() - 1)
				break;
			if (lightPath[i].direction.length() < 0.5f)
				break;

			vec3f scatterFactor = (lightPath[i].color * lightPath[i].getCosineTerm() / 
				(lightPath[i + 1].originProb * lightPath[i].directionProb));

			lightState.throughput *= scatterFactor;
			/*
			vec3f decay = lightPath[i].getRadianceDecay((lightPath[i].origin - lightPath[i + 1].origin).length());
			fprintf(fp , "l = %d , thr = (%.8f,%.8f,%.8f) , color = (%.8f,%.8f,%.8f)\ncosine = %.8f , dirPdf = %.8f , oPdf = %.8f\ndecay=(%.8f,%.8f,%.8f)\n" , 
				i , lightState.throughput.x , lightState.throughput.y , lightState.throughput.z ,
				lightPath[i].color[0] , lightPath[i].color[1] , lightPath[i].color[2] ,
				lightPath[i].getCosineTerm() , lightPath[i].directionProb , lightPath[i + 1].originProb ,
				decay.x , decay.y , decay.z);
			*/
			if (lightPath[i].directionSampleType == Ray::RANDOM && useWeight)
			{
				Real pdf = lightPath[i].directionProb;
				if (pdf < 1e-7f)
					break;

				Real weightFactor;

				Real volMergeScale = 1;
				Real originProb = 1.f / totArea;
				Real dirProb;
				//dirProb = INV_2_PI;
				dirProb = lightPath[i].getCosineTerm() / M_PI;

				//if (lightPath[i].insideObject && lightPath[i].contactObject)
				//	printf("!!!\n");
				if (lightPath[i].insideObject && !lightPath[i].contactObject && lightPath[i].insideObject->isVolumetric())
				{
					volMergeScale = 4.f / 3.f * mergeRadius;
					originProb = 1.f / totVol;
					dirProb = 0.25f / M_PI;
				}

				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &dirProb));

				if (_isnan(weightFactor) || abs(pdf) < 1e-6f || abs(mergeFactor(&volMergeScale , &originProb , &dirProb)) < 1e-6f)
				{
					printf("sample light path error, %.8f , %.8f\n" , connectFactor(pdf) , 
						mergeFactor(&volMergeScale , &originProb , &dirProb));
				}
				/*
				if (abs(volMergeScale - 1.f) < 1e-6)
					printf("surface %.8f\n" , weightFactor);
				else
					printf("volume %.8f %.8f\n" , weightFactor);
				*/
				lightState.throughput *= weightFactor;
			}
		}
	}

	lightPhotonNum = partialPhotonNum = partialSubPathList.size();
}

Ray IptTracer::genIntermediateSamples(Scene& scene)
{
	if (totVol > 1e-7f)
	{
		float volProb = 0.6f , surProb = 1 - volProb;
		if (totArea < 1e-7f)
		{
			volProb = 1.f; surProb = 0.f;
		}
		if (RandGenerator::genFloat() < volProb)
		{
			Ray ray = genVolumeSample();
			ray.originProb *= volProb;
			return ray;
		}
		else
		{
			Ray ray = genOtherSurfaceSample();
			ray.originProb *= surProb;
			return ray;
		}
	}
	else
	{
		Ray ray = genOtherSurfaceSample();
		return ray;
	}
}

void IptTracer::genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>& interPathList)
{
#pragma omp parallel for
	for(int p=0; p<interPathNum; p++)
	{
		if (!renderer->scene.usingGPU())
		{
			Ray interRay = genIntermediateSamples(renderer->scene);
			interPathList[p] = new Path;
			samplePath(*interPathList[p] , interRay);
		}

		Path& interPath = *interPathList[p];
		/*
		if (interPath.size() >= 3)
		{
			if (interPath[2].contactObject && interPath[2].contactObjectTriangleID >= interPath[2].contactObject->faceVertexIndexList.size())
			{
				printf("init contact triangle error, %d %d\n" , interPath[2].contactObjectTriangleID , interPath[2].contactObject->faceVertexIndexList.size());
			}
		}
		*/

		//fprintf(fp , "=================\n");
		partPathMergeIndex[p].clear();

		if (interPath.size() <= 1)
			continue;

		IptPathState interState;
		interState.originRay = &interPath[0];

		interState.throughput = interPath[0].color * interPath[0].getCosineTerm() / 
			(interPath[0].originProb * interPath[0].directionProb * interPath[1].originProb);
		interState.indirContrib = vec3f(0.f);

		//if (intensity(interState.throughput) > 30.f)
		//	continue;

		/*
		fprintf(fp , "====================\n");
		vec3f decay = interPath[0].getRadianceDecay((interPath[0].origin - interPath[1].origin).length());
		fprintf(fp , "l = 0 , thr = (%.8f,%.8f,%.8f) , color = (%.8f,%.8f,%.8f)\ncosine = %.8f , dirPdf = %.8f , oPdf = %.8f\ndecay=(%.8f,%.8f,%.8f)\n" , 
			interState.throughput.x , interState.throughput.y , interState.throughput.z ,
			interPath[0].color[0] , interPath[0].color[1] , interPath[0].color[2] ,
			interPath[0].getCosineTerm() , interPath[0].directionProb , interPath[1].originProb ,
			decay.x , decay.y , decay.z);
		*/

		for(unsigned i = 1; i < interPath.size(); i++)
		//for (unsigned i = 1; i < 2; i++)
		{
			Real dist = std::max((interPath[i].origin - interPath[i - 1].origin).length() , 1e-5f);
			interState.throughput *= interPath[i - 1].getRadianceDecay(dist);

			if(interPath[i].contactObject && interPath[i].contactObject->emissive())
				break;

			interState.pos = interPath[i].origin;
			interState.lastRay = &interPath[i - 1];
			interState.ray = &interPath[i];
			interState.pathLen = i;
			
			if(interPath[i].directionSampleType != Ray::DEFINITE &&
				(interPath[i].insideObject != NULL || interPath[i].contactObject != NULL) &&
				(interPath[i].origin != interPath[i - 1].origin))
				//(interPath[i].insideObject && !interPath[i].contactObject)) // only volume
			{
				//fprintf(fp , "path length = %d, dirContrib = (%.8f,%.8f,%.8f)\n" , 
				//	i , interState.dirContrib[0] , interState.dirContrib[1] , interState.dirContrib[2]);

				omp_set_lock(&cmdLock);
				partialSubPathList.push_back(interState);
				partPathMergeIndex[p].push_back(partialSubPathList.size() - 1);
				omp_unset_lock(&cmdLock);
			}

			if (i == interPath.size() - 1)
				break;
			if (interPath[i].direction.length() < 0.5f)
				break;

			/*
			if (interPath[i].contactObjectTriangleID >= interPath[i].contactObject->faceVertexIndexList.size())
				printf("contact triangle error, %d %d\n" , interPath[i].contactObjectTriangleID , interPath[i].contactObject->faceVertexIndexList.size());
			*/
			vec3f scatterFactor = (interPath[i].color * interPath[i].getCosineTerm() / 
				(interPath[i + 1].originProb * interPath[i].directionProb));

			interState.throughput *= scatterFactor;
			
			/*
			vec3f decay = interPath[i].getRadianceDecay((interPath[i].origin - interPath[i + 1].origin).length());
			fprintf(fp , "l = %d , thr = (%.8f,%.8f,%.8f) , color = (%.8f,%.8f,%.8f)\ncosine = %.8f , dirPdf = %.8f , oPdf = %.8f\ndecay=(%.8f,%.8f,%.8f)\n" , 
				i , interState.throughput.x , interState.throughput.y , interState.throughput.z ,
				interPath[i].color[0] , interPath[i].color[1] , interPath[i].color[2] ,
				interPath[i].getCosineTerm() , interPath[i].directionProb , interPath[i + 1].originProb ,
				decay.x , decay.y , decay.z);
			*/

			if (interPath[i].directionSampleType != Ray::DEFINITE && useWeight)
			{
				Real pdf = interPath[i].directionProb;
				if (pdf < 1e-7f)
					break;

				Real weightFactor;

				Real volMergeScale = 1.f;
				Real originProb = 1.f / totArea;
				Real dirProb;
				//dirProb = INV_2_PI;
				dirProb = interPath[i].getCosineTerm() / M_PI;

				//if (interPath[i].insideObject && interPath[i].contactObject)
				//	printf("!!!\n");
				if (interPath[i].insideObject && !interPath[i].contactObject && interPath[i].insideObject->isVolumetric())
				{
					volMergeScale = 4.f / 3.f * mergeRadius;
					originProb = 1.f / totVol;
					dirProb = 0.25f / M_PI;
				}
				
				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &dirProb));

				if (_isnan(weightFactor) || abs(pdf) < 1e-6f || abs(mergeFactor(&volMergeScale , &originProb , &dirProb)) < 1e-6f)
				{
					printf("sample inter path error, %.8f , %.8f\n" , connectFactor(pdf) , 
						mergeFactor(&volMergeScale , &originProb , &dirProb));
				}

				interState.throughput *= weightFactor;
			}
		}
	}

	partialPhotonNum = partialSubPathList.size();
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

	vector<vec3f> contribs(partialPhotonNum - lightPhotonNum);

	for (int i = 0; i < partialSubPathList.size(); i++)
	{
		partialSubPathList[i].index = i;
	}

	if (usePPM)
		return;

	// preprocess
	PointKDTree<IptPathState> lightTree(partialSubPathList);

	revIndex = new int[partialPhotonNum - lightPhotonNum];

	int flag = 0;

	for (int i = 0; i < interPathNum; i++)
	{
		if (partPathMergeIndex[i].size() == 0)
			continue;

		for (int j = 0; j < partPathMergeIndex[i].size(); j++)
		{
			int k = partPathMergeIndex[i][j] - lightPhotonNum;
			revIndex[k] = i;
		}
	}

#pragma omp parallel for
	for (int i = 0; i < interPathNum; i++)
	{
		if (partPathMergeIndex[i].size() == 0)
			continue;

		SearchQuery query;

		query.interState = &partialSubPathList[partPathMergeIndex[i][0]];

		lightTree.searchInRadius(0 , query.interState->originRay->origin , mergeRadius , query);

		partPathMergeIndex[i].clear();

		for (int j = 0; j < query.mergeIndex.size(); j++)
		{
			int k = query.mergeIndex[j];
			if (k < lightPhotonNum || revIndex[k - lightPhotonNum] != i)
				partPathMergeIndex[i].push_back(query.mergeIndex[j]);
		}
	}

	bool f = true;
	while (f)
	{
		f = false;
		// check cycle
		vis.resize(partialPhotonNum - lightPhotonNum);
		for (int i = 0; i < vis.size(); i++)
			vis[i] = 0;
		int totColor = 0;
		for (int st = lightPhotonNum; st < partialPhotonNum; st++)
		{
			if (vis[st - lightPhotonNum] > 0)
				continue;
			totColor++;
			while (!cycle.empty())
				cycle.pop();
			f |= dfs(st , totColor);
		}

		// eliminate cycle
		for (int i = 0; i < edgeToRemove.size(); i++)
		{
			int child = edgeToRemove[i].second;
			int pa = edgeToRemove[i].first;
			int k = revIndex[child - lightPhotonNum];
			int index = -1;
			for (int j = 0; j < partPathMergeIndex[k].size(); j++)
			{
				if (partPathMergeIndex[k][j] == pa)
					index = j;
			}
			if (index != -1)
				partPathMergeIndex[k].erase(partPathMergeIndex[k].begin() + index);
		}
		edgeToRemove.clear();
	}

	printf("check & eliminate cycles\n");
	
	f = true;
	int totMergeIter = 0;
	while (f)
	{
		++totMergeIter;
#pragma omp parallel for
		for (int i = lightPhotonNum; i < partialPhotonNum; i++)
		{
			mergePartialPaths(contribs , partialSubPathList[i]);
		}

		f = false;
		for (int i = lightPhotonNum; i < partialPhotonNum; i++)
		{
			if (!f && abs(intensity(partialSubPathList[i].indirContrib) - 
				intensity(contribs[i - lightPhotonNum])) > 1e-4f)
				f = true;
			partialSubPathList[i].indirContrib = contribs[i - lightPhotonNum];
		}

		//if (totMergeIter > mergeIterations) break;
	}

	printf("merge done... totMergeIter = %d... tracing eye paths...\n" , totMergeIter);

	vis.clear();
	for (int i = 0; i < partPathMergeIndex.size(); i++)
	{
		partPathMergeIndex[i].clear();
		partPathMergeIndex[i].shrink_to_fit();
	}
	partPathMergeIndex.clear();

	delete[] revIndex;
}

bool IptTracer::dfs(int cur , int col)
{
	if (partialSubPathList[cur].ray->contactObject)
		return false;
	if (vis[cur - lightPhotonNum] == col)
	{
		stack<int> tmp(cycle);
		while (!tmp.empty())
		{
			if (cur == tmp.top())
			{
				// mark cycle edges
				edgeToRemove.push_back(pair<int , int>(cur , cycle.top()));
				// end of mark

				// print cycle
				/*
				cycle.push(cur);
				fprintf(fp , "!!!!!!!!!!!!!!!! CYCLE !!!!!!!!!!!!!!!\n");
				while (!cycle.empty())
				{
					fprintf(fp , "%d " , cycle.top());
					int x = cycle.top();
					cycle.pop();
					int y;
					if (cycle.empty())
						break;
					else
						y = cycle.top();
					fprintf(fp , "%.8f " , (partialSubPathList[x].pos - 
						partialSubPathList[y].originRay->origin).length());
				}
				fprintf(fp , "\n");
				*/
				// end of print

				return true;
			}
			tmp.pop();
		}
		return false;
	}
	cycle.push(cur);
	vis[cur - lightPhotonNum] = col;
	int k = revIndex[cur - lightPhotonNum];
	bool isCycle = false;
	for (int j = 0; j < partPathMergeIndex[k].size(); j++)
	{
		int pa = partPathMergeIndex[k][j];
		if (pa < lightPhotonNum)
			continue;

		isCycle |= dfs(pa , col);
		if (isCycle)
			return true;
	}
	cycle.pop();
	return false;
}

vec3f IptTracer::colorByConnectingLights(Ray lastRay , Ray ray)
{
	Ray lightRay = genEmissiveSurfaceSample();
	lightRay.direction = (ray.origin - lightRay.origin);
	Real dist = lightRay.direction.length();
	Real dist2 = dist * dist;
	lightRay.direction.normalize();
	Ray outRay = ray;
	outRay.direction = -lightRay.direction;

	vec3f decayFactor = outRay.getRadianceDecay(dist);

	if(!testVisibility(outRay, lightRay))
		return vec3f(0.f);

	//outRay.direction = -cameraState.lastRay->direction;
	//vec3f bsdfFactor2 = lightRay.getBSDF(outRay);
	vec3f bsdfFactor = lastRay.getBSDF(outRay);

	if (y(bsdfFactor) < 1e-7f)
		return vec3f(0.f);

	Real cosAtLight = clampf(lightRay.getContactNormal().dot(lightRay.direction) , 0.f , 1.f);

	Real cosToLight = clampf(ray.getContactNormal().dot(-lightRay.direction) , 0.f , 1.f);

	if (cosAtLight < 1e-6f && cosToLight < 1e-6f)
		return vec3f(0.f);

	vec3f tmp = lightRay.color * cosAtLight * bsdfFactor * cosToLight
		/ (lightRay.originProb * dist2);

	Real p1 = lightRay.originProb;
	Real p2 = lightRay.originProb * (cosAtLight / M_PI) * cosToLight / dist2 * 
		(partialPathNum * mergeRadius * mergeRadius * M_PI);
	Real weightFactor = p1 / (p1 + p2);

	//fprintf(fp , "weight = %.8f , bsdfToLightPdf = %.8f , cosAtLight = %.8f ,\ntoLightOriginPdf = %.8f , originProb = %.8f , dist = %.8f\n" , 
	//	weightFactor , bsdfToLightPdf , cosAtLight , toLightOriginPdf , lightRay.originProb , dist);

	vec3f res = tmp * decayFactor;

	res *= weightFactor;
	
	//printf("sur weight = %.8f\n" , weightFactor);

	/*
	vec3f resx = camera.eliminateVignetting(res , cameraState.index) * pixelNum;
	if (cameraState.ray->contactObject)//(resx[0] + resx[1] + resx[2] >= 2)
	{
		fprintf(fp , "=====================\n");
		fprintf(fp , "cosAtLight = %.8f, cosToLight = %.8f, originPdf = %.8f, pdf = %.8f, weight=%.8f,\nres=(%.10f,%.10f,%.10f)\nbsdf=(%.10f,%.10f,%.10f)\n" , 
			cosAtLight , cosToLight , originProb , pdf , weightFactor , resx[0] , resx[1] , resx[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2]);
	}
	*/
	return res;
}

vec3f IptTracer::colorByRayMarching(Path& eyeMergePath , PointKDTree<IptPathState>& partialSubPaths)
{
	vec3f tr(1.f) , surfaceRes(0.f) , volumeRes(0.f);
	for (int i = 1; i < eyeMergePath.size(); i++)
	{
		// hit light
		if (eyeMergePath[i].contactObject && eyeMergePath[i].contactObject->emissive())
		{
			surfaceRes = eyeMergePath[i].color;
			break;
		}

		// volume
		Real dist = std::max((eyeMergePath[i - 1].origin - eyeMergePath[i].origin).length() , 1e-5f);
		if (eyeMergePath[i - 1].insideObject && eyeMergePath[i - 1].insideObject->isVolumetric())
		{
			if (eyeMergePath[i - 1].insideObject->isHomogeneous())
			{
				GatherQuery query(this);
				query.color = vec3f(0.f);
				query.constKernel = false;

				Ray volRay = eyeMergePath[i - 1];
				SceneVPMObject *vol = static_cast<SceneVPMObject*>(volRay.insideObject);
				Real stepSize = vol->stepSize;
				int N = dist / stepSize;
				if (N == 0)
					N++;

				Real step = dist / N;
				Real offset = step * RandGenerator::genFloat();
				float t = offset;
				tr *= vol->getRadianceDecay(volRay , offset);

				IptPathState cameraState;
				cameraState.throughput = vec3f(1.f);
				cameraState.lastRay = &volRay;

				Ray outRay = volRay;
				outRay.contactObject = NULL;

				for (int i = 0; i < N; i++)
				{
					query.color = vec3f(0.f);
					outRay.origin = volRay.origin + volRay.direction * t;
					cameraState.ray = &outRay;
					cameraState.pos = outRay.origin;
					query.cameraState = &cameraState;

					partialSubPaths.searchInRadius(0 , query.cameraState->pos , mergeRadius , query);

					outRay.origin -= volRay.direction * t;
					tr *= vol->getRadianceDecay(outRay , step);

					volumeRes += query.color * tr * step;
					t += step;
				}
			}
			else
			{
				GatherQuery query(this);
				query.color = vec3f(0.f);
				query.constKernel = false;

				Ray volRay = eyeMergePath[i - 1];
				HeterogeneousVolume *vol = static_cast<HeterogeneousVolume*>(volRay.insideObject);
				float stepSize = vol->getStepSize();
				int N = dist / stepSize;
				if (N == 0)
					N++;

				float step = dist / N;
				float offset = step * RandGenerator::genFloat();
				float t = offset;
				tr *= vol->getRadianceDecay(volRay , offset);

				IptPathState cameraState;
				cameraState.throughput = vec3f(1.f);
				cameraState.lastRay = &volRay;

				Ray outRay = volRay;
				outRay.contactObject = NULL;

				for(int i = 0; i < N; i++)
				{
					query.color = vec3f(0.f);
					outRay.origin = volRay.origin + volRay.direction * t;
					cameraState.ray = &outRay;
					cameraState.pos = outRay.origin;
					query.cameraState = &cameraState;

					partialSubPaths.searchInRadius(0 , query.cameraState->pos , mergeRadius , query);

					outRay.origin -= volRay.direction * t;
					tr *= vol->getRadianceDecay(outRay , step);

					volumeRes += query.color * tr * step;
					
					t += step;
				}
			}
		}
		else
		{
			if (eyeMergePath[i - 1].insideObject)
				tr *= eyeMergePath[i - 1].getRadianceDecay(dist);
		}

		// surface
		if (eyeMergePath[i].contactObject && eyeMergePath[i].directionSampleType == Ray::RANDOM)
		{
			if (eyeMergePath[i].contactObject->isVolumetric())
				continue;

			GatherQuery query(this);
			query.color = vec3f(0.f);
			query.constKernel = false;

			IptPathState cameraState;
			cameraState.throughput = vec3f(1.f);
			cameraState.lastRay = &eyeMergePath[i - 1];
			cameraState.ray = &eyeMergePath[i];
			cameraState.pos = cameraState.ray->origin;

			query.cameraState = &cameraState;

			partialSubPaths.searchInRadius(0 , query.cameraState->pos , mergeRadius , query);

			surfaceRes = query.color;

			if (!usePPM)
			{
				vec3f dirIllu = colorByConnectingLights(eyeMergePath[i - 1] , eyeMergePath[i]);
				surfaceRes += dirIllu;
			}

			break;
		}
	}

	vec3f color = tr * surfaceRes + volumeRes;
	return color;
}

void IptTracer::mergePartialPaths(vector<vec3f>& contribs , const IptPathState& interState)
{
	struct MergeQuery
	{
		vec3f color;
		IptTracer *tracer;
		const IptPathState* interState;

		MergeQuery(IptTracer* tracer) { this->tracer = tracer; }

		void process(const IptPathState& lightState)
		{
			Real dist = (lightState.ray->origin - interState->originRay->origin).length();
			Real volMergeScale = 1.f;
			if (interState->originRay->insideObject && !interState->originRay->contactObject && 
				interState->originRay->insideObject->isVolumetric())
				volMergeScale = 4.f / 3.f * tracer->mergeRadius;
			
			if (dist >= tracer->mergeRadius)
				printf("dist error\n");

			if (interState->originRay->insideObject && !interState->originRay->contactObject &&
				interState->originRay->insideObject->isVolumetric())
			{
				if (lightState.ray->insideObject == NULL || 
					lightState.ray->contactObject != NULL ||
					interState->originRay->insideObject != lightState.ray->insideObject ||
					!interState->originRay->insideObject->isVolumetric() ||
					!lightState.ray->insideObject->isVolumetric())
				{
					return;
				}
				
				volMergeScale = 4.0 / 3.0 * tracer->mergeRadius;
			}
			else if (interState->originRay->contactObject)
			{
				//if (interState->originRay->insideObject)
				//	printf("should not have insideObject\n");

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
			
			vec3f totContrib;
			if (lightState.index < tracer->lightPhotonNum)
				totContrib = lightState.throughput;
			else
				totContrib = lightState.indirContrib;

			if (y(totContrib) < 1e-7f)
				return;

			Ray outRay , inRay;
			vec3f bsdfFactor;
		
			outRay = *lightState.ray;
			outRay.direction = interState->originRay->direction;

			inRay = *lightState.lastRay;
			inRay.direction = interState->originRay->origin - lightState.lastRay->origin;
			inRay.direction.normalize();
			//vec3f bsdfFactor2 = lightState.lastRay->getBSDF(outRay);
			bsdfFactor = inRay.getBSDF(*interState->originRay);

			//if ((bsdfFactor - bsdfFactor2).length() > 1e-6f)
			//{
			//	printf("bsdf error, (%.8f,%.8f,%.8f),(%.8f,%.8f,%.8f)\n" , 
			//		bsdfFactor.x , bsdfFactor.y , bsdfFactor.z ,
			//		bsdfFactor2.x , bsdfFactor2.y , bsdfFactor2.z);
			//}

			if (y(bsdfFactor) < 1e-7f)
				return;

			vec3f tmp = totContrib * bsdfFactor * interState->throughput;

			Real lastPdf , weightFactor;
			//Real lastPdf2 = lightState.lastRay->getDirectionSampleProbDensity(outRay);
			lastPdf = inRay.getDirectionSampleProbDensity(*interState->originRay);

			//if (abs(lastPdf - lastPdf2) > 1e-6f)
			//{
			//	printf("pdf error, %.8f, %.8f\n" , lastPdf , lastPdf2);
			//}

			//if (lastPdf < 1e-7f)
			//	return;

			/*
			fprintf(fp , "================\n");
			fprintf(fp , "f1=(%.8f,%.8f,%.8f), f2=(%.8f,%.8f,%.8f), p1=%.8f, p2=%.8f\n" ,
				bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , bsdfFactor2[0] , bsdfFactor2[1] , bsdfFactor2[2] , 
				lastPdf , lastPdf2);
			*/
			
			weightFactor = tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb) /
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb));

			if (_isnan(weightFactor) || abs(lastPdf) < 1e-6f || 
				abs(tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb)) < 1e-6f)
			{
				printf("merge partial path error, %.8f , %.8f\n" , tracer->connectFactor(lastPdf) , 
					tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb));
			}

			vec3f res;
			res = tmp * (tracer->mergeKernel / volMergeScale);

			//printf("%.8f , %.8f\n" , interState->originRay->originProb , interState->originRay->directionProb);
			
			res *= weightFactor;

			color += res;
			/*
			vec3f resx = res;	
			if (resx.x > totContrib.x || resx.y > totContrib.y || resx.z > totContrib.z)
			{
				fprintf(fp , "----------------\n");
				fprintf(fp , "res = (%.8f,%.8f,%.8f), totContrib = (%.8f,%.8f,%.8f), \nbsdf = (%.8f,%.8f,%.8f), \ninterThr = (%.8f,%.8f,%.8f), weightFactor = %.8f\nPc = %.8f, Pm = %.8f, originProb = %.8f, dirProb = %.8f\n" ,
					resx[0] , resx[1] , resx[2] , totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , 
					interState->throughput[0] , interState->throughput[1] , interState->throughput[2] , weightFactor ,
					tracer->connectFactor(lastPdf) , tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb ,
					&interState->originRay->directionProb) , interState->originRay->originProb , interState->originRay->directionProb);
			}
			*/
		}
	};

	MergeQuery query(this);
	query.interState = &interState;
	query.color = vec3f(0, 0, 0);

	int pa = revIndex[interState.index - lightPhotonNum];
	/*
	if (partPathMergeIndex[pa].size() > 0)
		fprintf(fp , "======%d=======\n" , partPathMergeIndex[pa].size());
	*/

	for (int j = 0; j < partPathMergeIndex[pa].size(); j++)
	{
		int k = partPathMergeIndex[pa][j];
		query.process(partialSubPathList[k]);
	}

	contribs[interState.index - lightPhotonNum] = query.color;
}

void IptTracer::movePaths(omp_lock_t& cmdLock , vector<Path>& pathListGPU , vector<Path*>& pathList)
{
	for (int i = 0; i < pathListGPU.size(); i++)
	{
		pathList[i] = &pathListGPU[i];
	}
}

void IptTracer::sampleMergePath(Path &path, Ray &prevRay, uint depth) 
{
	path.push_back(prevRay);

	Ray terminateRay;
	terminateRay.origin = prevRay.origin;
	terminateRay.color = vec3f(0,0,0);
	terminateRay.direction = vec3f(0,0,0);
	terminateRay.directionSampleType = Ray::DEFINITE;
	terminateRay.insideObject =	terminateRay.contactObject = terminateRay.intersectObject = NULL;

	Ray nextRay;
	if (prevRay.insideObject && !prevRay.insideObject->isVolumetric())
	{
		nextRay = prevRay.insideObject->scatter(prevRay);
	}
	else 
	{
		if (prevRay.intersectObject)
		{
			if (prevRay.intersectObject->isVolumetric() && 
				prevRay.contactObject && prevRay.contactObject->isVolumetric())
			{
				prevRay.origin += prevRay.direction * prevRay.intersectDist;
				prevRay.intersectDist = 0;
			}
			nextRay = prevRay.intersectObject->scatter(prevRay);
		}
		else
		{
			path.push_back(terminateRay);
			return;
		}
	}

	if (nextRay.direction.length() < 0.5) 
	{
		path.push_back(nextRay);		
		return;
	}

	if (depth + 1 > maxDepth)
	{
		path.push_back(terminateRay);
		return;
	}
	
	terminateRay.origin = nextRay.origin;
	if (nextRay.contactObject && (nextRay.contactObject->emissive() || nextRay.directionSampleType == Ray::RANDOM))
	{
		path.push_back(nextRay);
		path.push_back(terminateRay);
		return;
	}
	
	NoSelfIntersectionCondition condition(&renderer->scene, nextRay);
	Scene::ObjSourceInformation info;
	float dist = renderer->scene.intersect(nextRay, info, &condition);
	if (dist < 0)
	{
		path.push_back(nextRay);
		path.push_back(terminateRay);
		return;
	}
	else
	{
		nextRay.intersectObject = renderer->scene.objects[info.objID];
		nextRay.intersectObjectTriangleID = info.triangleID;
		nextRay.intersectDist = dist;
	}
	sampleMergePath(path, nextRay, depth + 1);
}