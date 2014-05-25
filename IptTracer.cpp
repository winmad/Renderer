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
#pragma omp parallel for
			for (int i = 0; i < partialPhotonNum; i++)
			{
				IptPathState& subPath = partialSubPathList[i];

				int _x(0) , _y(0);
				vec3f color(0.f);
				color = colorByConnectingCamera(camera , subPath , _x , _y);
				
				if (y(color) > 0)
				{
					omp_set_lock(&pixelLocks[_y*camera.width + _x]);
					singleImageColors[_y*camera.width + _x] += color * 0.5f;
					omp_unset_lock(&pixelLocks[_y*camera.width + _x]);
				}
			}
			*/
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
			/*
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

			if (totVol > 1e-7f)
			{
				//countHashGrid.print(fp);
			}

			if (!usePPM)
			{
				vector<IptPathState> lightSubPathList(partialSubPathList);
#pragma omp parallel for
				for (int p = 0; p < interPathNum; p++)
					interRayList[p] = genIntermediateSamples(lightSubPathList , renderer->scene);
				interPathListGPU = samplePathList(interRayList);
				movePaths(cmdLock , interPathListGPU , interPathList);

				genIntermediatePaths(cmdLock , interPathList);
			}
			
			printf("partialPhotonNum = %d\n" , partialSubPathList.size());

			mergeKernel = 1.f / (M_PI * mergeRadius * 
				mergeRadius * (Real)partialPathNum);

			if (!usePPM)
				mergePartialPaths(cmdLock);

			PointKDTree<IptPathState> partialSubPaths(partialSubPathList);
			
			printf("merge done... tracing eye paths...\n");

#pragma omp parallel for
			for (int p = 0; p < cameraPathNum; p++)
				eyeRayList[p] = camera.generateRay(p);
			eyePathListGPU = samplePathList(eyeRayList);

			for(int p=0; p<cameraPathNum; p++)
			{
				Path eyePath = eyePathListGPU[p];

				if (eyePath.size() <= 1)
					continue;

				IptPathState cameraState;
				cameraState.isSpecularPath = 1;

				cameraState.throughput = vec3f(1.f) / (eyePath[0].originProb * eyePath[0].directionProb * eyePath[1].originProb);
	
				cameraState.index = eyePath.front().pixelID;

				//cameraState.accuProb = 1.f;

				vector<float> ratios , weights;
				//calcEyeProbRatios(eyePath , ratios);

				//fprintf(fp , "===================\n");

				int nonSpecLength = 0;
				vector<vec3f> mergeContribs;
				vec3f mergeContrib(0.f);

				mergeContribs.clear();
				weights.clear();
				weights.push_back(1.f);

				for(unsigned i = 1; i < eyePath.size(); i++)
				//for (unsigned i = 1; i < 3; i++)
				{
					vec3f colorHitLight(0.f) , colorDirIllu(0.f) , colorGlbIllu(0.f);

					Real dist = std::max((eyePath[i].origin - eyePath[i - 1].origin).length() , 1e-5f);
					cameraState.throughput *= eyePath[i - 1].getRadianceDecay(dist);

					if(eyePath[i].contactObject && eyePath[i].contactObject->emissive())
					{
						if (cameraState.isSpecularPath)
						{
							vec3f le = ((SceneEmissiveObject*)(eyePath[i].contactObject))->getColor();
							colorHitLight = le * cameraState.throughput;
							singleImageColors[cameraState.index] += colorHitLight;
						}	
						break;
					}

					if ((eyePath[i].insideObject == NULL && eyePath[i].contactObject == NULL) ||
						(eyePath[i].origin == eyePath[i - 1].origin))
						break;

					cameraState.pos = eyePath[i].origin;
					cameraState.lastRay = &eyePath[i - 1];
					cameraState.ray = &eyePath[i];

					Real eyeWeight = 1.f;
					//eyeWeight = calcEyePathWeight(eyePath , ratios , i);

					if(eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						//colorDirIllu = colorByConnectingLights(camera , cameraState);
						//singleImageColors[cameraState.index] += colorDirIllu;

						colorGlbIllu = colorByMergingPaths(cameraState , partialSubPaths);
						mergeContribs.push_back(colorDirIllu + colorGlbIllu);
						//weights.push_back(cameraState.accuProb);
					}

					if (eyePath[i].directionSampleType == Ray::RANDOM)
					{
						cameraState.isSpecularPath = 0;
						nonSpecLength++;

						//if (usePPM)
						if (nonSpecLength == 1)
							break; // PPM, eye path length is 1
					}

					if (i >= eyePath.size() - 1)
						break;
				
					vec3f bsdfFactor;
					bsdfFactor = eyePath[i].color;
				
					cameraState.throughput *= (bsdfFactor * eyePath[i].getCosineTerm() / 
						(eyePath[i + 1].originProb * eyePath[i].directionProb));

					//fprintf(fp , "l=%d, thr=(%.8f,%.8f,%.8f), bsdf=(%.8f,%.8f,%.8f), cos=%.8f, prob=%.8f\n" , 
					//	i , cameraState.throughput[0] , cameraState.throughput[1] , cameraState.throughput[2] ,
					//	bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , eyePath[i].getCosineTerm() , eyePath[i].directionProb);
					
					if (eyePath[i].directionSampleType != Ray::DEFINITE)
					{
						Real pdf = eyePath[i].directionProb;
		
						Real weightFactor;

						Real volMergeScale = 1;
						Real originProb = 1.f / totArea;
						//originProb = getOriginProb(countHashGrid , cameraState.pos , false);
						if (eyePath[i].insideObject && eyePath[i].contactObject == NULL)
						{
							volMergeScale = 4.0 / 3.0 * mergeRadius;
							originProb = 1.f / totVol;
							//originProb = getOriginProb(countHashGrid , cameraState.pos , true);
						}

						//originProb = lastAccuProb;
						
						weightFactor = connectFactor(pdf) /
							(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &INV_2_PI));

						//fprintf(fp , "w = %.8f\n" , weightFactor);

						//cameraState.throughput *= weightFactor;
						//weights.push_back(weightFactor);
					}
				}
				
				Real sumWeights = 0.f;
				for (int i = 0; i < weights.size(); i++)
					sumWeights += weights[i];
				for (int i = 0; i < mergeContribs.size(); i++)
				{
					singleImageColors[cameraState.index] += mergeContribs[i] * weights[i] / sumWeights;
					//if (nonSpecLength > 0)
					//	singleImageColors[cameraState.index] += mergeContribs[i] / (Real)nonSpecLength;		
				}	
			}
			*/
		}

		printf("done calculation, release memory\n");
		// add to image & release memory
		if(cmd == "exit")
			return pixelColors;

		//eliminateVignetting(singleImageColors);

		for(int i=0; i<pixelColors.size(); i++)
		{
			pixelColors[i] *= (Real)s / ((Real)s + 1.f);
			pixelColors[i] += singleImageColors[i] / ((Real)s + 1.f);// * camera.width * camera.height; 
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
		lightState.isSpecularPath = 1;
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

				if (_isnan(weightFactor))
				{
					printf("%.8f , %.8f\n" , connectFactor(pdf) , 
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
		float volProb = 1.f , surProb = 1 - volProb;
		if (totArea < 1e-7f)
		{
			volProb = 1.f; surProb = 0.f;
		}
		if (RandGenerator::genFloat() < volProb)
		{
			Ray ray = genVolumeSample();
			ray.originProb *= volProb;
			//printf("%.8f, %.8f\n" , ray.originProb , ray.directionProb);
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
		//printf("%.8f , %.8f\n" , ray.originProb , ray.directionProb);
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
		interState.isSpecularPath = 1;
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

			interState.isSpecularPath &= (interPath[i].directionSampleType == Ray::DEFINITE);
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

				if (_isnan(weightFactor))
				{
					printf("%.8f , %.8f\n" , connectFactor(pdf) , 
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

	// check
	/*
	for (int i = 0; i < interPathNum; i++)
	{
		//int id = (int)(RandGenerator::genFloat() * interPathNum);
		int id = i;
		if (partPathMergeIndex[id].size() == 0)
			continue;

		int cnt = 0;
		IptPathState interState = partialSubPathList[partPathMergeIndex[id][0]];

		//fprintf(fp , "===========*==========\n");
		for (int j = 0; j < partPathMergeIndex[id].size(); j++)
		{
			IptPathState state = partialSubPathList[partPathMergeIndex[id][j]];
			//fprintf(fp , "origin=(%.8f,%.8f,%.8f), pos=(%.8f,%.8f,%.8f), thr=(%.8f,%.8f,%.8f)\noriginProb = %.8f, dirProb = %.8f\n" ,
			//	state.originRay->origin.x , state.originRay->origin.y , state.originRay->origin.z ,
			//	state.pos.x , state.pos.y , state.pos.z ,
			//	state.throughput.x , state.throughput.y , state.throughput.z ,
			//	state.originRay->originProb , state.originRay->directionProb);

			if (partPathMergeIndex[id][j] < lightPhotonNum)
				printf("error partPathMergeIndex\n");
			if (interState.originRay != partialSubPathList[partPathMergeIndex[id][j]].originRay)
				printf("error originRay\n");
		}
	}
	*/

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

		//fprintf(fp , "===========\n");
		//fprintf(fp , "interPos = (%.8f,%.8f,%.8f)\n" , query.interState->pos[0] , 
		//	query.interState->pos[1] , query.interState->pos[2]);

		lightTree.searchInRadius(0 , query.interState->originRay->origin , mergeRadius , query);

		/*
		fprintf(fp , "=======mergeNum = %d========\n" , query.mergeIndex.size());
		for (int j = 0; j < query.mergeIndex.size(); j++)
		{
			fprintf(fp , "%d " , query.mergeIndex[j]);
		}
		fprintf(fp , "\n");
		*/

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
			int k = query.mergeIndex[j];
			if (k < lightPhotonNum || revIndex[k - lightPhotonNum] != i)
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

	// check
	/*
	for (int i = 0; i < partialPhotonNum - lightPhotonNum; i++)
	{
		//int id = (int)(RandGenerator::genFloat() * (partialPhotonNum - lightPhotonNum));
		int id = i;
		int k = revIndex[id];

		if (partPathMergeIndex[k].size() == 0)
			continue;

		int cnt = 0;
		IptPathState interState = partialSubPathList[lightPhotonNum + id];

		fprintf(fp , "=*******************=\n");
		for (int j = 0; j < partPathMergeIndex[k].size(); j++)
		{
			vec3f dis = interState.originRay->origin - partialSubPathList[partPathMergeIndex[k][j]].pos;
			fprintf(fp , "%.8f " , dis.length());
		}
		fprintf(fp , "\n");
	
		fprintf(fp , "=====================\n");
		fprintf(fp , "%d\n" , partPathMergeIndex[k].size());
		for (int j = 0; j < partPathMergeIndex[k].size(); j++)
			fprintf(fp , "%d " , partPathMergeIndex[k][j]);
		fprintf(fp , "\n");

		for (int j = 0; j < partialPhotonNum; j++)
		{
			vec3f dis = interState.originRay->origin - partialSubPathList[j].pos;

			if (dis.length() <= mergeRadius)
			{
				fprintf(fp , "%d " , j);
				cnt++;
			}
		}
		fprintf(fp , "\n%d\n" , cnt);
		if (cnt != partPathMergeIndex[k].size())
			printf("error count\n");
	}
	*/

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
		dfs(st , totColor);
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

	printf("check & eliminate cycles\n");
	
	bool f = true;
	int totMergeIter = 0;
	//for (int mergeIter = 0; mergeIter < mergeIterations; mergeIter++)
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
				intensity(contribs[i - lightPhotonNum])) > 1e-6f)
				f = true;
			partialSubPathList[i].indirContrib = contribs[i - lightPhotonNum];
		}
	}

	printf("merge done... totMergeIter = %d... tracing eye paths...\n" , totMergeIter);

	edgeToRemove.clear();
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

vec3f IptTracer::colorByConnectingCamera(const Camera& camera, const IptPathState& lightState , int& _x , int& _y)
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
	Real originProb = 1.f / totArea;
	//originProb = getOriginProb(countHashGrid , lightState.ray->origin , false);
	if (lightState.ray->insideObject && lightState.ray->contactObject == NULL)
	{
		volMergeScale = 4.0 / 3.0 * mergeRadius;
		originProb = 1.f / totVol;
		//originProb = getOriginProb(countHashGrid , lightState.ray->origin , true);
	}

	Real imagePointToCameraDist = camera.sightDist / cosAtCamera;
	Real imageToSolidAngleFactor = imagePointToCameraDist *
		imagePointToCameraDist / cosAtCamera;
	Real imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / distEye2;

	Real cameraPdfArea = imageToSurfaceFactor * 1.f; // pixel area is 1
	
	Real surfaceToImageFactor = 1.f / imageToSurfaceFactor;
	
	Real bsdfDirPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);

	vec3f totContrib = lightState.throughput + lightState.indirContrib;

	//---- still buggy, fix me ----
	vec3f color = (totContrib * bsdfFactor * decayFactor) //* cosAtCamera * cosToCamera / distEye2;
		 / (surfaceToImageFactor * lightPathNum);
	color *= powf(cosAtCamera , 4.f) / pixelNum;
	//-----------------------------

	Ray inRay;
	inRay.origin = camera.position;
	inRay.direction = -dirToCamera;

	if (!testVisibility(inRay , outRay))
		return vec3f(0.f);

	Real pdf = bsdfDirPdf;

	Real weightFactor = connectFactor(pdf) / 
		(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &INV_2_PI));

	//fprintf(fp , "weight = %.8f\n" , weightFactor);

	if (useWeight)
		color *= weightFactor;

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

vec3f IptTracer::colorByConnectingLights(const Camera& camera, const IptPathState& cameraState)
{
	Ray lightRay = genEmissiveSurfaceSample();
	lightRay.direction = (cameraState.pos - lightRay.origin);
	Real dist = lightRay.direction.length();
	Real dist2 = dist * dist;
	lightRay.direction.normalize();
	Ray outRay = *cameraState.ray;
	outRay.direction = -lightRay.direction;
	Ray inRay = *cameraState.lastRay;

	vec3f decayFactor = outRay.getRadianceDecay(dist);

	if(!testVisibility(outRay, lightRay))
		return vec3f(0.f);

	//outRay.direction = -cameraState.lastRay->direction;
	//vec3f bsdfFactor2 = lightRay.getBSDF(outRay);
	vec3f bsdfFactor = cameraState.lastRay->getBSDF(outRay);

	if (y(bsdfFactor) < 1e-7f)
		return vec3f(0.f);

	Real cosAtLight = min(max(0.f , lightRay.getContactNormal().dot(lightRay.direction)) , 1.f);

	Real cosToLight = 1;
	if (cameraState.ray->contactObject)
	{
		cosToLight = min(max(0.f , cameraState.ray->getContactNormal().dot(-lightRay.direction)) , 1.f);
	}

	vec3f tmp = lightRay.color * cosAtLight * bsdfFactor * cosToLight * cameraState.throughput
		/ (lightRay.originProb * dist2);

	Real bsdfToLightPdf = inRay.getDirectionSampleProbDensity(outRay);

	outRay.direction = -cameraState.lastRay->direction;
	//Real lightOriginPdf = lightRay.getOriginSampleProbDensity(outRay);
	Real lightOriginPdf = lightRay.originProb;
	Real pdf = lightRay.getDirectionSampleProbDensity(outRay);

	Real volMergeScale = 1;
	Real originProb = 1.f / totArea;
	//originProb = getOriginProb(countHashGrid , cameraState.ray->origin , false);
	if (cameraState.ray->insideObject && !cameraState.ray->contactObject)
	{
		volMergeScale = 4.0 / 3.0 * mergeRadius;
		originProb = 1.f / totVol;
		//originProb = getOriginProb(countHashGrid , cameraState.ray->origin , true);
	}

	//fprintf(fp , "p1=%.8f, p2=%.8f\n" , bsdfToLightPdf , pdf);

	//originProb = lightRay.directionProb * lightOriginPdf * cosToLight / (dist2);
	//originProb = pdf * lightOriginPdf * cosToLight / (dist2);

	Real weightFactor = connectFactor(pdf) /
		(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &INV_2_PI));

	//fprintf(fp , "weight = %.8f , bsdfToLightPdf = %.8f , cosAtLight = %.8f ,\ntoLightOriginPdf = %.8f , originProb = %.8f , dist = %.8f\n" , 
	//	weightFactor , bsdfToLightPdf , cosAtLight , toLightOriginPdf , lightRay.originProb , dist);

	vec3f res = tmp * decayFactor;

	if (useWeight)
		res *= weightFactor;
	
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

vec3f IptTracer::colorByMergingVolume(IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths , vec3f& tr)
{
	vec3f volumeRes(0.f);
	Real dist = std::max((cameraState.lastRay->origin - cameraState.ray->origin).length() , 1e-6f);
	if (cameraState.lastRay->insideObject && cameraState.lastRay->insideObject->isVolumetric())
	{
		if (cameraState.lastRay->insideObject->isHomogeneous())
		{
			GatherQuery query(this);
			query.cameraState = &cameraState;
			query.color = vec3f(0.f);
			query.constKernel = false;

			Ray volRay = *(cameraState.lastRay);
			SceneVPMObject *vol = static_cast<SceneVPMObject*>(volRay.insideObject);
			Real stepSize = vol->stepSize;
			int N = dist / stepSize;
			if (N == 0)
				N++;

			//printf("%d\n" , N);

			Real step = dist / N;
			Real offset = step * RandGenerator::genFloat();
			float t = offset;
			tr *= vol->getRadianceDecay(volRay , offset);

			Ray outRay = volRay;
			outRay.contactObject = NULL;
			query.cameraState->ray = &outRay;

			for (int i = 0; i < N; i++)
			{
				query.color = vec3f(0.f);
				query.cameraState->pos = volRay.origin + volRay.direction * t;
				
				outRay.origin = query.cameraState->pos;
				query.cameraState->ray->origin = query.cameraState->pos;
				query.cameraState->ray->direction = query.cameraState->lastRay->direction;

				partialSubPaths.searchInRadius(0 , query.cameraState->pos , mergeRadius , query);

				tr *= vol->getRadianceDecay(outRay , step);

				volumeRes += query.color * tr * step;
				t += step;
			}
		}
	}
	return volumeRes;
}

vec3f IptTracer::colorByMergingSurface(IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths)
{
	vec3f surfaceRes(0.f);
	if (cameraState.ray->contactObject && cameraState.ray->directionSampleType == Ray::RANDOM)
	{
		if (cameraState.ray->contactObject->isVolumetric())
			return surfaceRes;

		GatherQuery query(this);
		query.cameraState = &cameraState;
		query.color = vec3f(0.f);
		query.constKernel = false;

		partialSubPaths.searchInRadius(0 , query.cameraState->pos , mergeRadius , query);

		surfaceRes = query.color;
	}
	return surfaceRes;
}

vec3f IptTracer::colorByMergingPaths(IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths)
{	
	vec3f tr(1.f) , surfaceRes(0.f) , volumeRes(0.f);

	IptPathState tmpCameraState = cameraState;
	volumeRes = colorByMergingVolume(tmpCameraState , partialSubPaths , tr);

	tmpCameraState = cameraState;
	surfaceRes = colorByMergingSurface(tmpCameraState , partialSubPaths);

	cameraState.throughput *= tr;
	return (surfaceRes * tr + volumeRes);
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
			
			/*
			Real mergeKernel = tracer->kernel(dist * dist , tracer->mergeRadius * tracer->mergeRadius);
			if (abs(volMergeScale - 1.f) > 1e-6f)
				mergeKernel /= (float)tracer->partialPathNum * tracer->mergeRadius * tracer->mergeRadius * tracer->mergeRadius;
			else 
				mergeKernel /= (float)tracer->partialPathNum * tracer->mergeRadius * tracer->mergeRadius;
			*/

			weightFactor = tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb) /
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb));

			//Real mergeFactor = 1.f / mergeKernel * interState->originRay->originProb * interState->originRay->directionProb;
			//weightFactor = mergeFactor / (lastPdf + mergeFactor);

			if (_isnan(weightFactor))
			{
				printf("%.8f , %.8f\n" , tracer->connectFactor(lastPdf) , 
					tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb , &interState->originRay->directionProb));
			}

			vec3f res;
			res = tmp * (tracer->mergeKernel / volMergeScale);
			
			//res = tmp * mergeKernel;

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

Ray IptTracer::link(const Path& path, int i, int j)
{
	const Ray& src = path[clamp(i, 0, path.size()-1)];
	const Ray& dst = path[clamp(j, 0, path.size()-1)];
	Ray ray = src;
	ray.direction = dst.origin - src.origin;
	ray.direction.normalize();
	ray.insideObject = dst.insideObject;
	return ray;
}

void IptTracer::calcEyePathProbs(Path& eyePath , vector<double>& probDir , vector<double>& probRev)
{
	int N = eyePath.size() - 1;
	for (int i = N; i > 0; i--)
	{
		if (abs((eyePath[i].origin - eyePath[i - 1].origin).length()) < 1e-6f)
		{
			N = i;
			break;
		}
	}

	probDir.clear();
	probRev.clear();
	probDir.resize(N);
	probRev.resize(N);

	probDir.front() = eyePath.front().originProb;

	for(int i=0; i<N-1; i++)
	{
		probDir[i+1] = probDir[i] * link(eyePath, i+1, i).getCosineTerm();
		double dist = max2((eyePath[i+1].origin - eyePath[i].origin).length(), EPSILON);
		probDir[i+1] /= dist * dist;

		if(eyePath[i].directionSampleType == Ray::RANDOM)
		{
			if(i>0)
			{
				probDir[i+1] *= eyePath[i].directionProb;
			}
			else
			{
				probDir[i+1] *= 1.f;
			}
			probDir[i+1] *= eyePath[i+1].originProb;
		}

	}
	/*
	if (eyePath[N].insideObject && !eyePath[N].contactObject)
		probRev.back() = getOriginProb(countHashGrid , eyePath[N].origin , true);
	else
		probRev.back() = 1.f / totArea;
	*/
	probRev.back() = 1.f;

	for(int i = N-1; i>0; i--)
	{
		probRev[i-1] = probRev[i] * eyePath[i-1].getCosineTerm();
		double dist = max2((eyePath[i-1].origin - eyePath[i].origin).length(), EPSILON);
		probRev[i-1] /= dist * dist;

		if(eyePath[i].directionSampleType == Ray::RANDOM)
		{
			if(i < N)
			{
				probRev[i-1] *= link(eyePath, i+1, i).getDirectionSampleProbDensity(link(eyePath, i, i-1));
			}
			else
			{
				if (eyePath[N].insideObject && !eyePath[N].contactObject)
					probRev[i-1] *= 0.25f / M_PI;
				else
					probRev[i-1] *= 0.5f / M_PI;
			}
			probRev[i-1] *= link(eyePath, i, i-1).getOriginSampleProbDensity(link(eyePath, i-1, i-2));
		}
	}
	/*
	fprintf(fp , "======================\n");
	for (int i = 0; i < N; i++)
		fprintf(fp , "%.8f " , probDir[i]);
	fprintf(fp , "\n");
	for (int i = 0; i < N; i++)
		fprintf(fp , "%.8f " , probRev[i]);
	fprintf(fp , "\n");
	*/
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