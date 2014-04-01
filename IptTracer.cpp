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
	totVol = renderer->scene.getTotalVolume();
	printf("scene: totArea = %.8f, totVol = %.8f\n" , totArea , totVol);

	for(int i=0; i<pixelLocks.size(); i++)
	{
		omp_init_lock(&pixelLocks[i]);
	}

	omp_lock_t cmdLock;
	omp_init_lock(&cmdLock);

	Real r0 = mergeRadius;

	countHashGrid.init(renderer->scene);

	for(unsigned s=0; s<spp; s++)
	{
		if (!renderer->scene.usingGPU())
		{
			partPathMergeIndex.resize(lightPathNum + interPathNum);

			partialSubPathList.clear();

			mergeRadius = r0 * (pow(s+1, (alpha-1)/3.f));
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

			countHashGrid.addPhotons(partialSubPathList , 0 , lightPhotonNum);

			/*
			fprintf(fp , "============ one iter ============\n");
			for (int i = 0; i < countHashGrid.weights.size(); i++)
			{
				if (countHashGrid.weights[i] <= 0)
					continue;
				fprintf(fp , "i = %d, weight = %.8f\n" , i , countHashGrid.weights[i]);
			}
			*/

			genIntermediatePaths(cmdLock , interPathList);
	
			mergePartialPaths(cmdLock);

			countHashGrid.addPhotons(partialSubPathList , lightPhotonNum , partialPhotonNum);

			printf("partialPathNum = %d\n" , partialSubPathList.size());

			int st = 0;
			int ed = lightPhotonNum;
			int surNum = 0 , volNum = 0;
#pragma omp parallel for
			for (int i = st; i < ed; i++)
			{
				IptPathState& subPath = partialSubPathList[i];

				subPath.dirContrib = vec3f(0.f);

				int _x(0) , _y(0);
				vec3f color(0.f);
				//color = colorByConnectingCamera(camera , subPath , _x , _y);
				
				if (y(color) > 0)
				{
					omp_set_lock(&pixelLocks[_y*camera.width + _x]);
					singleImageColors[_y*camera.width + _x] += color;
					omp_unset_lock(&pixelLocks[_y*camera.width + _x]);
				}
				/*
				if (s == 0)
				{
					if (subPath.ray->insideObject && !subPath.ray->contactObject)
					{
						fprintf(fp , "volume, dirContrib=(%.4f,%.4f,%.4f), indirContrib=(%.4f,%.4f,%.4f), throughput=(%.4f,%.4f,%.4f)\n" ,
							subPath.dirContrib[0] , subPath.dirContrib[1] , subPath.dirContrib[2] ,
							subPath.indirContrib[0] , subPath.indirContrib[1] , subPath.indirContrib[2] ,
							subPath.throughput[0] , subPath.throughput[1] , subPath.throughput[2]);
					}
					else
					{
						fprintf(fp , "surface, dirContrib=(%.4f,%.4f,%.4f), indirContrib=(%.4f,%.4f,%.4f), throughput=(%.4f,%.4f,%.4f)\n" ,
							subPath.dirContrib[0] , subPath.dirContrib[1] , subPath.dirContrib[2] ,
							subPath.indirContrib[0] , subPath.indirContrib[1] , subPath.indirContrib[2] ,
							subPath.throughput[0] , subPath.throughput[1] , subPath.throughput[2]);
					}
				}
				*/
			}
			//printf("sur = %d, vol = %d\n" , surNum , volNum);

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

				cameraState.accuProb = initialProb / (eyePath[0].directionProb);
				//cameraState.accuProb = initialProb;

				//omp_set_lock(&cmdLock);
				int T = eyePath.size();
				for (int i = 1; i < T; i++)
				{
					if ((eyePath[i].contactObject && eyePath[i].contactObject->emissive()) ||
						(eyePath[i].insideObject == NULL && eyePath[i].contactObject == NULL) ||
						(eyePath[i].origin == eyePath[i - 1].origin))
					{
						T = i;
						break;
					}
				}
				vector<float> ratios;
				for (int i = 1; i < T; i++)
				{
					float r , dist , deno = 1.f , nume = 1.f;
					Ray inRay , outRay , tmpRay;
					vec3f dir;
					if (T < 3)
						break;
					if (i == 1)
					{
						inRay = eyePath[i + 1];
						inRay.direction = eyePath[i].origin - eyePath[i + 1].origin;
						inRay.direction.normalize();
						outRay = eyePath[i];
						outRay.direction = eyePath[i - 1].origin - eyePath[i].origin;
						dist = outRay.direction.length();
						outRay.direction.normalize();
						tmpRay = eyePath[i - 1];
						if (outRay.directionSampleType == Ray::DEFINITE)
							nume = eyePath[i].directionProb;
						else
							nume = inRay.getDirectionSampleProbDensity(outRay);
						//nume *= outRay.getOriginSampleProbDensity(tmpRay);
						nume *= eyePath[i].originProb;
						nume *= eyePath[i - 1].getCosineTerm() / (dist * dist);
						deno = 1.f;
						
						/*
						if (nume > 1e10f || abs(nume) < 1e-10f);
						{
							fprintf(fp , "=======nume error=======\n");
							fprintf(fp , "nume = %.8f, %.8f, %.8f , %.8f , %.8f\n" , nume , inRay.getDirectionSampleProbDensity(outRay) , 
								outRay.getOriginSampleProbDensity(tmpRay) , eyePath[i - 1].getCosineTerm() , dist);
						}
						*/
					}
					else if (i == T - 1)
					{
						// uniform origin prob
						if (eyePath[i].insideObject && !eyePath[i].contactObject)
							nume = 1.f / totVol;
						else
							nume = 1.f / totArea;
						dir = eyePath[i - 1].origin - eyePath[i].origin;
						dist = dir.length();
						dir.normalize();

						deno = (eyePath[i - 1].directionProb * eyePath[i].originProb
							 / (dist * dist));
						if (eyePath[i].getContactNormal().length() > 0)
							deno *= abs(eyePath[i].getContactNormal().dot(dir));
						/*
						if (nume > 1e6 || abs(nume) < 1e-6);
						{
							fprintf(fp , "=======nume error=======\n");
							fprintf(fp , "nume = %.8f\n" , nume);
						}

						if (deno > 1e6 || abs(deno) < 1e-6);
						{
							fprintf(fp , "=======deno error=======\n");
							fprintf(fp , "deno = %.8f, %.8f, %.8f , %.8f , %.8f\n" , deno , eyePath[i - 1].directionProb , eyePath[i].originProb , 
								abs(eyePath[i].getContactNormal().dot(dir)) , dist);
						}
						*/
					}
					else
					{
						if (i == T - 2)
						{
							outRay = eyePath[i + 1];
							outRay.direction = eyePath[i].origin - eyePath[i + 1].origin;
							dist = outRay.direction.length();
							outRay.direction.normalize();
							tmpRay = eyePath[i];
							nume = INV_2_PI * outRay.getOriginSampleProbDensity(tmpRay)
								* eyePath[i].getCosineTerm() / (dist * dist);				
						}
						else
						{
							inRay = eyePath[i + 2];
							inRay.direction = eyePath[i + 1].origin - eyePath[i + 2].origin;
							inRay.direction.normalize();
							outRay = eyePath[i + 1];
							outRay.direction = eyePath[i].origin - eyePath[i + 1].origin;
							dist = outRay.direction.length();
							outRay.direction.normalize();
							tmpRay = eyePath[i];
							if (eyePath[i + 1].directionSampleType == Ray::DEFINITE)
								nume = eyePath[i + 1].directionProb;
							else 
								nume = inRay.getDirectionSampleProbDensity(outRay);
							nume *= outRay.getOriginSampleProbDensity(tmpRay)
								* eyePath[i].getCosineTerm() / (dist * dist);
						}

						dir = eyePath[i - 1].origin - eyePath[i].origin;
						dist = dir.length();
						dir.normalize();
						deno = (eyePath[i - 1].directionProb * eyePath[i].originProb 
							 / (dist * dist));
						if (eyePath[i].getContactNormal().length() > 0)
							deno *= abs(eyePath[i].getContactNormal().dot(dir));
					}
					r = nume / deno;
					ratios.push_back(r);
				}

				/*
				fprintf(fp , "=========eye path=========\n");
				for (int i = 0; i < ratios.size(); i++)
				{
					fprintf(fp , "ratio = %.8f, weight = %.8f\n" , ratios[i] , 
						calcEyePathWeight(eyePath , ratios , i + 1));
				}
				*/

				//omp_unset_lock(&cmdLock);

				for(unsigned i=1; i<eyePath.size(); i++)
				{
					vec3f colorHitLight(0.f) , colorDirIllu(0.f) , colorGlbIllu(0.f);

					Real dist = std::max((eyePath[i].origin - eyePath[i - 1].origin).length() , 1e-5f);
					cameraState.throughput *= eyePath[i - 1].getRadianceDecay(dist);

					Real lastAccuProb = cameraState.accuProb;
					cameraState.accuProb *= eyePath[i - 1].directionProb * eyePath[i].originProb;
					if (std::abs(eyePath[i].getCosineTerm()) > 1e-5f)
						cameraState.accuProb /= std::abs(eyePath[i - 1].getCosineTerm());
					
					if(eyePath[i].contactObject && eyePath[i].contactObject->emissive())
					{
						if (cameraState.isSpecularPath)
						{
							vec3f le = ((SceneEmissiveObject*)(eyePath[i].contactObject))->getColor();
							colorHitLight = le * cameraState.throughput;
							omp_set_lock(&pixelLocks[cameraState.index]);
							//singleImageColors[cameraState.index] += colorHitLight;
							omp_unset_lock(&pixelLocks[cameraState.index]);
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
						colorGlbIllu = colorByMergingPaths(cameraState , partialSubPaths);
					}

					omp_set_lock(&pixelLocks[cameraState.index]);
					singleImageColors[cameraState.index] += colorDirIllu * eyeWeight;
					singleImageColors[cameraState.index] += colorGlbIllu * eyeWeight;
					omp_unset_lock(&pixelLocks[cameraState.index]);

					if (eyePath[i].contactObject != NULL && 
						(eyePath[i].directionSampleType == Ray::RANDOM || eyePath[i].contactObject->canMerge))
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
						Real originProb = 1.f / totArea;
						if (eyePath[i].insideObject && eyePath[i].contactObject == NULL)
						{
							volMergeScale = 4.0 / 3.0 * mergeRadius;
							originProb = 1.f / totVol;
						}

						//originProb = lastAccuProb;
						
						weightFactor = connectFactor(pdf) /
							(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &pdf));

						//fprintf(fp , "w = %.8f\n" , weightFactor);

						//cameraState.throughput *= weightFactor;
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

			if (clock() / 1000 >= lastTime)
			{
				showCurrentResult(pixelColors , &lastTime);
				lastTime += timeInterval;
			}
			else
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
			Real dist = std::max((lightPath[i].origin - lightPath[i - 1].origin).length() , 1e-5f);
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
			/*
			if (lightState.isSpecularPath)
			{
				fprintf(fp , "==========light path===========\n");
				fprintf(fp , "bsdf = (%.8f,%.8f,%.8f), pdf = %.8f, cos = %.8f, dirContrib = (%.8f,%.8f,%.8f), thr = (%.8f,%.8f,%.8f)\n" , 
					lightPath[i].color[0] , lightPath[i].color[1] , lightPath[i].color[2] , lightPath[i].directionProb ,
					lightPath[i].getCosineTerm() , lightState.dirContrib[0] , lightState.dirContrib[1] , lightState.dirContrib[2] ,
					lightState.throughput[0] , lightState.throughput[1] , lightState.throughput[2]);
			}
			*/
			if (lightPath[i].directionSampleType != Ray::DEFINITE)
			{
				Real weightFactor;

				Real volMergeScale = 1;
				Real originProb = 1.f / totArea;
				if (lightPath[i].insideObject && lightPath[i].contactObject == NULL)
				{
					volMergeScale = 4.0 / 3.0 * mergeRadius;
					originProb = 1.f / totVol;
				}

				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &pdf));

				lightState.throughput *= weightFactor;
				lightState.dirContrib *= weightFactor;
			}
		}
	}

	lightPhotonNum = partialSubPathList.size();
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
		float volScale = 0.f;
		if (lightSubPathList[i].ray->insideObject && !lightSubPathList[i].ray->contactObject)
			volScale = 1.f;
		weights[i + 1] = weights[i] + intensity * intensity * volScale;
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
			Real dist = std::max((interPath[i].origin - interPath[i - 1].origin).length() , 1e-5f);
			interState.throughput *= interPath[i - 1].getRadianceDecay(dist);

			if(interPath[i].contactObject && interPath[i].contactObject->emissive())
				break;

			interState.pos = interPath[i].origin;
			interState.lastRay = &interPath[i - 1];
			interState.ray = &interPath[i];
			
			if(interPath[i].directionSampleType != Ray::DEFINITE  &&
				(interPath[i].insideObject != NULL || interPath[i].contactObject != NULL) &&
				(interState.pos != interState.originRay->origin) &&
				(interPath[i].insideObject && !interPath[i].contactObject)) // only volume
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

			if (interPath[i].directionSampleType != Ray::DEFINITE)
			{
				Real weightFactor;

				Real volMergeScale = 1;
				Real originProb = 1.f / totArea;
				if (interPath[i].insideObject && interPath[i].contactObject == NULL)
				{
					volMergeScale = 4.0 / 3.0 * mergeRadius;
					originProb = 1.f / totVol;
				}

				//originProb = getOriginProb(countHashGrid , interPath[i].origin);
				//fprintf(fp , "originProb = %.8f\n" , originProb);
				
				weightFactor = connectFactor(pdf) /
					(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &pdf));
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

		/*
		fprintf(fp , "mergeNum = %d\n" , query.mergeIndex.size());
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
	for (;;) {
	float randWeight = RandGenerator::genFloat();
	int pathId = (lower_bound(weights.begin() , weights.end() , randWeight) - weights.begin()) - 1;
	pathId = clamp(pathId , 0 , partialSubPathList.size() - 1);
	IptPathState& lightState = partialSubPathList[pathId];
	while (lightState.ray == NULL ||
		(lightState.ray->insideObject && !lightState.ray->insideObject->canMerge) ||
		(!lightState.ray->insideObject && lightState.ray->contactObject && !lightState.ray->contactObject->canMerge))
	{
		randWeight = RandGenerator::genFloat();
		pathId = (lower_bound(weights.begin() , weights.end() , randWeight) - weights.begin()) - 1;
		pathId = clamp(pathId , 0 , partialSubPathList.size() - 1);
		lightState = partialSubPathList[pathId];
	}
	
	Ray ray;
	ray.originSampleType = Ray::SampleType::RANDOM;
	ray.directionSampleType = Ray::SampleType::RANDOM;

	ray.origin = lightState.pos;
	ray.insideObject = lightState.ray->insideObject;
	ray.contactObject = lightState.ray->contactObject;
	ray.contactObjectTriangleID = lightState.ray->contactObjectTriangleID;

	RandGenerator rng;
	Ray inRay , outRay(ray);
	inRay = *lightState.lastRay;

	vec3f o = lightState.pos;
	vec3f dir(0.f);

	//string msg("");

	if (lightState.ray->contactObject)
	{
		outRay = lightState.ray->contactObject->scatter(inRay , false);
		dir = outRay.direction;
		//msg = "surface";
	}
	else if (lightState.ray->insideObject)
	{
		outRay = lightState.ray->insideObject->scatter(inRay , false);
		dir = outRay.direction;
		//msg = "volume";
	}
	
	if (dir.length() < 1e-7f)
	{
		dir = rng.genSphericalDirection();
		outRay.directionProb = 0.25f / M_PI;
		//msg = "random";
	}
	dir.normalize();
	ray.origin = o + (dir * 1e-4 * rng.genFloat());

	ray.direction = dir;

	ray.color = vec3f(1.0);

	//ray.originProb = (weights[pathId + 1] - weights[pathId]);
	//ray.originProb = 1.f / totVol;
	ray.originProb = getOriginProb(countHashGrid , ray.origin);

	ray.directionProb = outRay.directionProb;

	vec3f bsdfFactor = inRay.getBSDF(outRay);
	if (y(bsdfFactor) < 1e-5)
		continue;

	/*
	fprintf(fp , "===============\n");
	fprintf(fp , "lightState: pos=(%.3f,%.3f,%.3f), c=(%.6f,%.6f,%.6f), bsdf=(%.6f,%.6f,%.6f)\n" , lightState.pos[0] , lightState.pos[1] ,
		lightState.pos[2] , lightState.dirContrib[0] , lightState.dirContrib[1] , lightState.dirContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2]);
	fprintf(fp , "------- %s -------\n" , msg.c_str());
	fprintf(fp , "interState: pos=(%.3f,%.3f,%.3f), dir=(%.3f,%.3f,%.3f), dirProb = %.6f, %.6f\n" , ray.origin[0] , ray.origin[1] , ray.origin[2] ,
		ray.direction[0] , ray.direction[1] , ray.direction[2] , ray.directionProb , inRay.getDirectionSampleProbDensity(outRay));
	*/
	
	ray.current_tid = scene.getContactTreeTid(ray);
	Scene::ObjSourceInformation osi;
	NoSelfIntersectionCondition condition(&scene , ray);
	Real dist = scene.intersect(ray, osi, &condition);

	if (dist > 0)
	{
		ray.intersectDist = dist;
		ray.intersectObject = scene.objects[osi.objID];
		ray.intersectObjectTriangleID = osi.triangleID;
	}
	if (ray.intersectObject == NULL)
		continue;

	return ray;
	}
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

	Real originProb = 1.f / totArea;
	Real volMergeScale = 1;
	if (lightState.ray->insideObject && lightState.ray->contactObject == NULL)
	{
		volMergeScale = 4.0 / 3.0 * mergeRadius;
		originProb = 1.f / totVol;
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
	//-----------------------------

	Ray inRay;
	inRay.origin = camera.position;
	inRay.direction = -dirToCamera;

	if (!testVisibility(inRay , outRay))
		return vec3f(0.f);

	Real pdf = bsdfDirPdf;
	//originProb = getOriginProb(countHashGrid , outRay.origin);

	Real weightFactor = connectFactor(pdf) / 
		(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &pdf));

	Real lightWeight = cameraPdfArea / pixelNum;
	lightWeight = 1.f / (1.f + lightWeight);

	//fprintf(fp , "weight = %.8f\n" , weightFactor);

	color *= lightWeight * weightFactor;
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
	//vec3f bsdfFactor = lightRay.getBSDF(outRay);
	vec3f bsdfFactor = cameraState.lastRay->getBSDF(outRay);
	if (y(bsdfFactor) < 1e-7f)
		return vec3f(0.f);

	Real cosAtLight = min(max(0.f , lightRay.getContactNormal().dot(lightRay.direction)) , 1.f);

	Real cosToLight = 1;
	if (cameraState.ray->contactObject != NULL)
	{
		cosToLight = min(max(0.f , cameraState.ray->getContactNormal().dot(-lightRay.direction)) , 1.f);
	}

	Real volMergeScale = 1;
	Real originProb = 1.f / totArea;
	if (cameraState.ray->insideObject && cameraState.ray->contactObject == NULL)
	{
		volMergeScale = 4.0 / 3.0 * mergeRadius;
		originProb = 1.f / totVol;
	}

	vec3f tmp = lightRay.color * cosAtLight * bsdfFactor * cosToLight * cameraState.throughput
		/ (lightRay.originProb * dist2);

	Real bsdfToLightPdf = inRay.getDirectionSampleProbDensity(outRay);

	outRay.direction = -cameraState.lastRay->direction;
	Real toLightOriginPdf = lightRay.getOriginSampleProbDensity(outRay);
	Real pdf = lightRay.getDirectionSampleProbDensity(outRay);

	originProb = cameraState.accuProb; // seems better with this line
	//originProb = getOriginProb(countHashGrid , outRay.origin);

	Real weightFactor = connectFactor(pdf) /
		(connectFactor(pdf) + mergeFactor(&volMergeScale , &originProb , &pdf));

	Real eyeWeight = bsdfToLightPdf / (lightRay.originProb / M_PI * cosToLight);

	eyeWeight = 1.f / (1.f + eyeWeight);

	/*
	fprintf(fp , "weight = %.8f , bsdfToLightPdf = %.8f , cosAtLight = %.8f ,\ntoLightOriginPdf = %.8f , originProb = %.8f , dist = %.8f\n" , 
		weightFactor , bsdfToLightPdf , cosAtLight , toLightOriginPdf , lightRay.originProb , dist);
	*/

	vec3f res = tmp * decayFactor * eyeWeight * weightFactor;
	
	/*
	vec3f resx = camera.eliminateVignetting(res , cameraState.index) * lightPathNum;
	if (resx[0] + resx[1] + resx[2] >= 2)
	{
		fprintf(fp , "=====================\n");
		fprintf(fp , "cosAtLight = %.8f, cosToLight = %.8f, originPdf = %.8f, pdf = %.8f, weight=%.8f,\nres=(%.10f,%.10f,%.10f)\nbsdf=(%.10f,%.10f,%.10f)\n" , 
			cosAtLight , cosToLight , originProb , pdf , weightFactor , resx[0] , resx[1] , resx[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2]);
	}
	*/
	return res;
}

vec3f IptTracer::colorByMergingPaths(const IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths)
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
			Real originProb = 1.f / tracer->totArea;
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
				originProb = 1.f / tracer->totVol;
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
			vec3f tmp = totContrib * bsdfFactor * cameraState->throughput; 

			Real lastPdf , weightFactor;
			lastPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);

			//originProb = cameraState->accuProb;
			//originProb = tracer->getOriginProb(tracer->countHashGrid , outRay.origin);

			weightFactor = tracer->mergeFactor(&volMergeScale , &originProb , &lastPdf) / 
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &originProb , &lastPdf));

			//fprintf(fp , "weight = %.8f , cameraAccuProb = %.8f , hashOriginProb = %.8f\n" , weightFactor , cameraState->accuProb , originProb);

			mergeNum++;
			vec3f res = tmp * (weightFactor * tracer->mergeKernel / volMergeScale);
			color += res;
			/*
			vec3f resx = tracer->renderer->camera.eliminateVignetting(res , cameraState->index) *
				tracer->pixelNum;	
			if (y(res) > 0)
			{
				fprintf(fp , "=====================\n");
				if (volMergeScale == 1)
					fprintf(fp , "surface\n");
				else 
					fprintf(fp , "volume\n");
				fprintf(fp , "res = (%.8f,%.8f,%.8f) \ntotContrib = (%.8f,%.8f,%.8f), bsdf = (%.8f,%.8f,%.8f),\n cameraThr = (%.8f,%.8f,%.8f) \nweightFactor = %.8f, originProb = %.8f, lastpdf = %.8f\n" ,
					resx[0] , resx[1] , resx[2] ,
					totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , 
					cameraState->throughput[0] , cameraState->throughput[1] , cameraState->throughput[2] , 
					weightFactor , originProb , lastPdf);
			}
			*/
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

			Real cosTerm = 1.f;
			if (lightState.ray->contactObject)
				cosTerm = lightState.ray->getContactNormal().dot(outRay.direction);

			vec3f tmp = totContrib * bsdfFactor * interState->throughput * cosTerm;

			Real lastPdf , weightFactor;
			lastPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);
			
			weightFactor = tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb ,
				&interState->originRay->directionProb) /
				(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &interState->originRay->originProb ,
				&interState->originRay->directionProb));
			
			//weightFactor = tracer->mergeFactor(&volMergeScale , NULL , &interState->originRay->directionProb) /
			//	(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , NULL , &interState->originRay->directionProb));

			color += tmp * (weightFactor * tracer->mergeKernel / volMergeScale);
			/*
			vec3f resx = tmp * weightFactor;	
			if (resx[0] + resx[1] + resx[2] > 0)
			{
				fprintf(fp , "=====================\n");
				fprintf(fp , "res = (%.8f,%.8f,%.8f), totContrib = (%.8f,%.8f,%.8f), \nbsdf = (%.8f,%.8f,%.8f), \ninterThr = (%.8f,%.8f,%.8f), weightFactor = %.8f, cosine = %.8f\n" ,
					resx[0] , resx[1] , resx[2] , totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , 
					interState->throughput[0] , interState->throughput[1] , interState->throughput[2] , weightFactor , cosTerm);
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