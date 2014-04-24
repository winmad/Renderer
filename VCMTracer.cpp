#include "StdAfx.h"
#include "VCMTracer.h"

static FILE* fp = fopen("debug_vcm.txt" , "w");

vector<vec3f> VCMTracer::renderPixels(const Camera& camera)
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

	float r0 = mergeRadius;

	for(unsigned s=0; s<spp; s++)
	{
		if(!renderer->scene.usingGPU())
		{

			mergeRadius = r0 * sqrt(powf(s+1, alpha-1));

			vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));

			string cmd;

			unsigned t = clock();

			vector<Path*> lightPathList(pixelColors.size(), NULL);

			vector<LightPathPoint> lppList;

#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Ray lightRay = genEmissiveSurfaceSample();
				lightPathList[p] = new Path;
				Path &lightPath = *lightPathList[p];
				samplePath(lightPath, lightRay);
				/*
				vec3f dirContrib(0.f);

				dirContrib = lightPath[0].color * lightPath[0].getCosineTerm() / 
					(lightPath[0].originProb * lightPath[0].directionProb *
					lightPath[1].originProb);
				
				if (s == 0)
					fprintf(fp , "=============\n");
				*/
				for(unsigned i=1; i<lightPath.size(); i++)
				{
					if(lightPath[i].contactObject && lightPath[i].contactObject->emissive())
						break;
					/*
					Real dist = std::max((lightPath[i].origin - lightPath[i - 1].origin).length() , 1e-5f);
					vec3f decayFactor = lightPath[i - 1].getRadianceDecay(dist);
					dirContrib *= decayFactor;
					
					if(lightPath[i].directionSampleType == Ray::RANDOM && s == 0 &&
						lightPath[i].insideObject && !lightPath[i].contactObject)
					{
						fprintf(fp , "path length = %d, dirContrib = (%.8f,%.8f,%.8f)\n" , 
							i , dirContrib[0] , dirContrib[1] , dirContrib[2]);
					}
					
					if (i + 1 < lightPath.size())
						dirContrib *= (lightPath[i].color * lightPath[i].getCosineTerm()) /
							(lightPath[i].directionProb * lightPath[i + 1].originProb);
					*/
					if(lightPath[i].directionSampleType == Ray::DEFINITE)
						continue;
					LightPathPoint lpp;
					lpp.pos = lightPath[i].origin;
					lpp.path = &lightPath;
					lpp.index = i;
					omp_set_lock(&cmdLock);
					lppList.push_back(lpp);
					omp_unset_lock(&cmdLock);
				}
			}

			PointKDTree<LightPathPoint> tree(lppList);

#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Path eyePath;

				const Path &lightPath = *lightPathList[p];

				samplePath(eyePath, camera.generateRay(p));

				colorByMergingPaths(singleImageColors, eyePath, tree);

				colorByConnectingPaths(pixelLocks, renderer->camera, singleImageColors, eyePath, lightPath);

			}

			if(cmd == "exit")
				return pixelColors;

			eliminateVignetting(singleImageColors);

			for(int i=0; i<pixelColors.size(); i++)
			{
				pixelColors[i] *= s / float(s + 1);
				pixelColors[i] += singleImageColors[i] / (s + 1)*camera.width*camera.height;
				delete lightPathList[i];
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
		else
		{
			/*
			mergeRadius = r0 * sqrt(powf(s+1, alpha-1));

			vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));

			string cmd;

			unsigned t = clock();

			vector<Path> eyePathList, lightPathList;
			vector<Ray> eyeRayList(pixelColors.size());
			vector<Ray> lightRayList(pixelColors.size());
#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				eyeRayList[p] = camera.generateRay(p);
				lightRayList[p] = genEmissiveSurfaceSample();
			}
			eyePathList = samplePathList(eyeRayList);
			lightPathList = samplePathList(lightRayList);
			vector<vector<unsigned>> visibilityList = testPathListVisibility(eyePathList, lightPathList);
			
			vector<LightPathPoint> lppList;

#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Path &lightPath = lightPathList[p];
				
				for(unsigned i=1; i<lightPath.size(); i++)
				{
					if(lightPath[i].contactObject && lightPath[i].contactObject->emissive())
						break;
					if(lightPath[i].directionSampleType == Ray::DEFINITE)
						continue;
					LightPathPoint lpp;
					lpp.pos = lightPath[i].origin;
					lpp.path = &lightPath;
					lpp.index = i;
					omp_set_lock(&cmdLock);
					lppList.push_back(lpp);
					omp_unset_lock(&cmdLock);
				}
			}

			PointKDTree<LightPathPoint> tree(lppList);

#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Path &eyePath = eyePathList[p];
				Path ep = eyePath;

				const Path &lightPath = lightPathList[p];

				colorByMergingPaths(singleImageColors, ep, tree);

				colorByConnectingPaths(pixelLocks, renderer->camera, singleImageColors, eyePath, lightPath, &visibilityList[p]);

			}

			if(cmd == "exit")
				return pixelColors;

			eliminateVignetting(singleImageColors);

			for(int i=0; i<pixelColors.size(); i++)
			{
				pixelColors[i] *= s / float(s + 1);
				pixelColors[i] += singleImageColors[i] / (s + 1)*camera.width*camera.height;
			}
			printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, (clock()-t_start)/1000);
			showCurrentResult(pixelColors);
			*/
		}
	}
	return pixelColors;
}

bool VCMTracer::mustUsePT(const Path& connectedPath)
{
	for(unsigned i=1; i<connectedPath.size()-1; i++)
		if(connectedPath[i].directionSampleType == Ray::RANDOM)
			return false;
	return true;
}

vec4<float> VCMTracer::connectColorProb(const Path& connectedPath, int connectIndex, bool merged)
{
	vec3f color(1, 1, 1);
	float prob = 1;

	float connectDist;

	for(int i=0; i<connectedPath.size(); i++)
	{
		color *= connectedPath[i].color;

		float dist = EPSILON;

		if(i <= connectIndex)
		{
			dist = max2((connectedPath[i+1].origin - connectedPath[i].origin).length(), EPSILON);
			color *= connectedPath[i].getRadianceDecay(dist);
		}
		else if(i>connectIndex+1)
		{
			dist = max2((connectedPath[i-1].origin - connectedPath[i].origin).length(), EPSILON);
			color *= connectedPath[i].getRadianceDecay(dist);
		}
		
		if(i==connectIndex && i<connectedPath.size()-1)
		{
			color *= connectedPath[i].getCosineTerm() * connectedPath[i+1].getCosineTerm() / (dist*dist);
			connectDist = dist;
		}

		if(i!=connectIndex && i!=connectIndex+1)
			color *= connectedPath[i].getCosineTerm();

		prob *= connectedPath[i].directionProb * connectedPath[i].originProb;
		if(merged && i == connectIndex)
		{
			prob *= M_PI * mergeRadius * mergeRadius;
			prob *= link(connectedPath, connectIndex, i+1, i).getCosineTerm();
			Ray prevRay = link(connectedPath, connectIndex, i, i+1);
			Ray nextRay = link(connectedPath, connectIndex, i+1, i+2);
			prob *= prevRay.getOriginSampleProbDensity(nextRay);
			if(connectedPath[i+1].insideObject && !connectedPath[i+1].contactObject)
			{
				prob *= 4.0/3.0 * mergeRadius;
			}
			prob /= dist * dist;
		}

		if(i > 0 && i <= connectIndex) // correct sample density difference in interpolating normal.
		{
			if (fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal())) != 0 &&
				fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal(true))) != 0)
			{
				color *= fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal())) /
					fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal(true)));
			}
			/*
			if (fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal())) == 0 ||
				fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal(true))) == 0)
			{
				fprintf(fp , "cos = %.8f , cos' = %.8f\n" , fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal())) ,
					fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal(true))));
				printf("error\n");
			}
			*/
		}
	}

	return vec4<float>(color, prob);
}

Ray VCMTracer::link(const Path& connectedPath, int connectIndex, int i, int j)
{
	const Ray& src = connectedPath[clamp(i, 0, connectedPath.size()-1)];
	const Ray& dst = connectedPath[clamp(j, 0, connectedPath.size()-1)];
	Ray ray = src;
	if((j>i && i<=connectIndex) || (j<i && i>connectIndex))
		return ray;
	ray.direction = dst.origin - src.origin;
	ray.direction.normalize();
	ray.insideObject = dst.insideObject;
	return ray;
}

float VCMTracer::connectMergeWeight(const Path& connectedPath, int connectIndex, bool merged, float expTerm)
{
	double sumExpProb = 0;

	vector<double> p_forward(connectedPath.size(), 1);
	vector<double> p_backward(connectedPath.size(), 1);
	vector<double> dist(connectedPath.size(), 0);

	p_forward.front() = connectedPath.front().getOriginSampleProbDensity(connectedPath.front());

	for(int i=0; i<connectedPath.size()-1; i++)
	{
		p_forward[i+1] = p_forward[i] * link(connectedPath, connectIndex, i+1, i).getCosineTerm();
		dist[i] = max2((connectedPath[i+1].origin - connectedPath[i].origin).length(), EPSILON);

		p_forward[i+1] /= dist[i]*dist[i];
		if(connectedPath[i].directionSampleType == Ray::RANDOM)
		{
			if(i>0)
			{
				p_forward[i+1] *= link(connectedPath, connectIndex, i-1, i).getDirectionSampleProbDensity(link(connectedPath, connectIndex, i, i+1));
			}
			else
			{
				Ray ray = link(connectedPath, connectIndex, i, i+1);
				p_forward[i+1] *= ray.getDirectionSampleProbDensity(ray);
			}
			p_forward[i+1] *= link(connectedPath, connectIndex, i, i+1).getOriginSampleProbDensity(link(connectedPath, connectIndex, i+1, i+2));
		}
	}

	p_backward.back() = connectedPath.back().getOriginSampleProbDensity(connectedPath.back());

	for(int i = connectedPath.size()-1; i>0; i--)
	{
		p_backward[i-1] = p_backward[i] * link(connectedPath, connectIndex, i-1, i).getCosineTerm() * connectedPath[i].originProb;
		p_backward[i-1] /= dist[i-1]*dist[i-1];
		if(connectedPath[i].directionSampleType == Ray::RANDOM)
		{
			if(i < connectedPath.size()-1)
			{
				p_backward[i-1] *= link(connectedPath, connectIndex, i+1, i).getDirectionSampleProbDensity(link(connectedPath, connectIndex, i, i-1));
			}
			else
			{
				Ray ray = link(connectedPath, connectIndex, i, i-1);
				p_backward[i-1] *= ray.getDirectionSampleProbDensity(ray);
			}
			p_backward[i-1] *= link(connectedPath, connectIndex, i, i-1).getOriginSampleProbDensity(link(connectedPath, connectIndex, i-1, i-2));
		}
	}

	double selfProb = connectIndex == -1 ? p_backward.front() : p_forward[connectIndex] * p_backward[connectIndex+1];

	for(int i=0; i<connectedPath.size()-1; i++)
	{
		if(connectedPath[i].directionSampleType == Ray::RANDOM && connectedPath[i+1].directionSampleType == Ray::RANDOM)
		{
			double p = p_forward[i] * p_backward[i+1];
			if(i == connectedPath.size() - 2)
				p *= renderer->camera.width * renderer->camera.height;
			sumExpProb += pow(p, double(expTerm));
		}
		if(i < connectedPath.size()-2 && connectedPath[i+1].directionSampleType == Ray::RANDOM)
		{
			double p = p_forward[i] * p_backward[i+1];
			if(connectedPath[i].directionSampleType == Ray::RANDOM)
			{
				if(i > 0)
				{
					p *= link(connectedPath, connectIndex, i-1, i).getDirectionSampleProbDensity(link(connectedPath, connectIndex, i, i+1));
				}
				else
				{
					Ray ray = link(connectedPath, connectIndex, i, i+1);
					p *= ray.getDirectionSampleProbDensity(ray);
				}
			}
			p *= link(connectedPath, connectIndex, i, i+1).getOriginSampleProbDensity(link(connectedPath, connectIndex, i+1, i+2));
			p *= M_PI * mergeRadius * mergeRadius;
			p *= link(connectedPath, connectIndex, i+1, i).getCosineTerm();
			p /= dist[i] * dist[i];
			if(connectedPath[i+1].insideObject && !connectedPath[i+1].contactObject)
			{
				p *= 4.0/3.0 * mergeRadius;
			}

			if(merged && connectIndex == i)
				selfProb = p;

			p *= renderer->camera.width * renderer->camera.height;
			sumExpProb += pow(p, double(expTerm));
		}
	}

	if(usePT || mustUsePT(connectedPath))
		sumExpProb += pow(p_backward.front(), double(expTerm));
	if(!(pow(selfProb, double(expTerm)) / sumExpProb >= 0))
	{
		return 0;
	}
	return pow(selfProb, double(expTerm)) / sumExpProb;
}

void VCMTracer::colorByConnectingPaths(vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const Path& eyePath, const Path& lightPath, vector<unsigned>* visibilityList)
{
	unsigned maxEyePathLen = eyePath.size();
	unsigned maxLightPathLen = lightPath.size();
	unsigned maxWholePathLen = maxEyePathLen + maxLightPathLen;

	for(int wholePathLen=2; wholePathLen<=maxWholePathLen; wholePathLen++)
	{

		Path wholePath(wholePathLen);

		int lplStart = (wholePathLen - int(maxEyePathLen)) > 0 ? (wholePathLen - maxEyePathLen) : 0;
		int lplEnd = wholePathLen - 1 < maxLightPathLen ? wholePathLen - 1 : maxLightPathLen;


		for(int lightPathLen = lplStart; lightPathLen <= lplEnd; lightPathLen++)
		{
			int eyePathLen = wholePathLen - lightPathLen;

			for(int li=0; li < lightPathLen; li++)
				wholePath[li] = lightPath[li];

			for(int ei=0; ei < wholePathLen - lightPathLen; ei++)
				wholePath[wholePathLen - 1 - ei] = eyePath[ei];

			int lightConnectID = lightPathLen-1;
			int eyeConnectID = lightPathLen;

			Ray& lightRay = wholePath[lightConnectID];
			Ray& eyeRay = wholePath[eyeConnectID];

			if(lightPathLen==0 && (!wholePath.front().contactObject || !wholePath.front().contactObject->emissive()))
				continue;

			if(lightPathLen==0 && !(usePT || mustUsePT(wholePath)))
				continue;

			if(eyePathLen > 0 && lightPathLen > 0)
			{
				if(lightPathLen > 1 && lightRay.contactObject && lightRay.contactObject->emissive())
					continue;
				if(eyeRay.contactObject && eyeRay.contactObject->emissive())
					continue;

				if(!connectRays(wholePath, lightConnectID))
					continue;

				if(!visibilityList)
				{
					if(!testVisibility(eyeRay, lightRay))
						continue;
				}
				else
				{
					unsigned visID = (eyePathLen-1)*maxLightPathLen + (lightPathLen-1);
					bool visible = ((*visibilityList)[visID/32] >> (visID%32)) & 1;
					if(!visible)
						continue;
				}
			}

			vec4<float> color_prob = connectColorProb(wholePath, lightConnectID);
				
			if(!(color_prob.w > 0))
				continue;

			if(!(vec3f(color_prob).length())>0)
				continue;

			float weight = connectMergeWeight(wholePath, lightConnectID, false);

			vec3f color = vec3f(color_prob) / color_prob.w * weight;

			if(!(color.length() > 0))
				continue;

			Ray &camRay = wholePath[wholePath.size()-1];

			if(eyePathLen > 1)// && eyePath[1].directionSampleType == Ray::DEFINITE)
			{
				omp_set_lock(&pixelLocks[camRay.pixelID]);

				//fprintf(fp , "w = %.8f\n" , weight);

				colors[camRay.pixelID] += color;
				omp_unset_lock(&pixelLocks[camRay.pixelID]);
			}
			else
			{
				/*
				if (eyePathLen <= 1)
				{
					fprintf(fp , "===========\n");
					for (int i = 0; i < wholePath.size(); i++)
					{
						fprintf(fp , "c=(%.8f,%.8f,%.8f), p=%.8f\n" , wholePath[i].color[0] , wholePath[i].color[1] ,
							wholePath[i].color[2] , wholePath[i].directionProb);
					}
				}
				*/
				vec2<float> pCoord = camera.transToPixel(wholePath[wholePath.size()-2].origin);
				int x = pCoord.x;
				int y = pCoord.y;
				if(x >= 0 && x < camera.width && y >= 0 && y < camera.height)
				{
					//fprintf(fp , "%.8f\n" , weight);
					omp_set_lock(&pixelLocks[y*camera.width + x]);
					colors[y*camera.width + x] += color;
					omp_unset_lock(&pixelLocks[y*camera.width + x]);
				}
			}
		}
	}
}

void VCMTracer::colorByMergingPaths(vector<vec3f>& colors, const Path& eyePath, PointKDTree<LightPathPoint>& tree)
{

	struct Query
	{
		vec3f color;
		VCMTracer *tracer;
		const Path* eyePath;
		unsigned eyePathLen;
		Path wholePath;

		int cnt;

		Query(VCMTracer* tracer) { this->tracer = tracer; }

		void process(const LightPathPoint& lpp)
		{
			Ray cameraRay = (*eyePath)[eyePathLen - 1];
			Ray lightRay = (*lpp.path)[lpp.index];
			if (cameraRay.insideObject && cameraRay.contactObject == NULL)
			{
				if (lightRay.insideObject == NULL ||
					lightRay.contactObject != NULL ||
					cameraRay.insideObject != lightRay.insideObject ||
					!cameraRay.insideObject->canMerge ||
					!lightRay.insideObject->canMerge)
				{
					return;
				}
			}
			else if (cameraRay.contactObject)
			{
				if (lightRay.contactObject != cameraRay.contactObject ||
					!lightRay.contactObject->canMerge ||
					!cameraRay.contactObject->canMerge)
				{
					return;
				}
			}
			else
			{
				return;
			}

			wholePath.assign(lpp.path->begin(), lpp.path->begin()+lpp.index);
			for(unsigned i=0; i<eyePathLen; i++)
			{
				wholePath.push_back((*eyePath)[eyePathLen-i-1]);
			}
			tracer->connectRays(wholePath, lpp.index-1, true);
			vec4<float> color_prob = tracer->connectColorProb(wholePath, lpp.index-1, true);
			if(vec3f(color_prob).length()==0 || color_prob.w==0)
				return;
			float weight = tracer->connectMergeWeight(wholePath, lpp.index-1, true);
			vec3f res = vec3f(color_prob) / color_prob.w * weight;

			color += res;
			
			//fprintf(fp , "res = (%.7f,%.7f,%.7f), c = (%.7f,%.7f,%.7f), prob = %.7f, weight = %.7f\n" ,
			//	res[0] , res[1] , res[2] , color_prob[0] , color_prob[1] , color_prob[2] , color_prob.w , weight);
			cnt++;
			
		}
	};

	Query query(this);
	query.eyePath = &eyePath;

	for(unsigned ei=1; ei<eyePath.size(); ei++)
	{
		if(eyePath[ei].contactObject && eyePath[ei].contactObject->emissive())
			break;
		if(eyePath[ei].directionSampleType == Ray::DEFINITE)
			continue;
		query.color = vec3f(0, 0, 0);
		query.eyePathLen = ei + 1;
		query.cnt = 0;

		tree.searchInRadius(0, eyePath[ei].origin, mergeRadius, query);
		//if (query.cnt > 0)
		//	fprintf(fp , "===================\n");

		colors[eyePath.front().pixelID] += query.color;
	}
}



