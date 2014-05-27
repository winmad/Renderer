#include "StdAfx.h"
#include "BidirectionalPathTracer.h"

static FILE* fp = fopen("debug_bpt.txt" , "w");

vector<vec3f> BidirectionalPathTracer::renderPixels(const Camera& camera)
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

	vector<vec3f> last_singleImageColors(pixelColors.size(), vec3f(0, 0, 0));
	vector<vec3f> colorWeights(pixelColors.size(), vec3f(0, 0, 0));

	for(unsigned s=0; s<spp; s++)
	{
		vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));

		string cmd;

		unsigned t = clock();

		if(!renderer->scene.usingGPU())
		{
#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Path eyePath;
				Path lightPath;

				Ray lightRay = genEmissiveSurfaceSample();

				samplePath(eyePath, camera.generateRay(p));

				samplePath(lightPath, lightRay);

				if(p == pathPixelID)
					showPath = lightPath;

				colorByConnectingPaths(pixelLocks, renderer->camera, singleImageColors, eyePath, lightPath);
			}
		}
		else
		{
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

#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Path& eyePath = eyePathList[p];
				Path& lightPath = lightPathList[p];

				colorByConnectingPaths(pixelLocks, renderer->camera, singleImageColors, eyePath, lightPath, &visibilityList[p]);

			}
		}

		if(cmd == "exit")
			return pixelColors;

		eliminateVignetting(singleImageColors);

		for(int i=0; i<pixelColors.size(); i++)
		{
			pixelColors[i] *= s / float(s + 1);
			pixelColors[i] += singleImageColors[i] / (s + 1);//*camera.width*camera.height;
		}

		printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, (clock()-t_start)/1000);

		//if (clock() / 1000 >= lastTime)
		if (s % outputIter == 0)
		{
			unsigned nowTime = clock() / 1000;
			showCurrentResult(pixelColors , &nowTime , &s);
			//showCurrentResult(pixelColors , &lastTime , &s);
			//lastTime += timeInterval;
		}
		else
			showCurrentResult(pixelColors);
	}
	return pixelColors;
}

bool BidirectionalPathTracer::mustUsePT(const Path& connectedPath)
{
	for(unsigned i=0; i<connectedPath.size()-1; i++)
		// changed
		if(connectedPath[i].directionSampleType == connectedPath[i+1].directionSampleType && connectedPath[i].directionSampleType == Ray::RANDOM)
			return false;
	return true;
}

vec4f BidirectionalPathTracer::connectColorProb(const Path& connectedPath, int connectIndex)
{
	vec3f color(1, 1, 1);
	float prob = 1;

	//fprintf(fp , "=====================\n");

	for(int i=0; i<connectedPath.size(); i++)
	{
		float dist;

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
		}

		if(i!=connectIndex && i!=connectIndex+1)
			color *= connectedPath[i].getCosineTerm();
		/*
		if(i > 0 && i <= connectIndex) // correct sample density difference in interpolating normal.
		{
			if (fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal())) != 0 &&
				fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal(true))) != 0)
			{
				color *= fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal())) /
					fabs(connectedPath[i-1].direction.dot(connectedPath[i].getContactNormal(true)));
			}
		}
		*/
		prob *= connectedPath[i].directionProb * connectedPath[i].originProb;

		color *= connectedPath[i].color;
		/*
		fprintf(fp , "len=%d, bsdf=(%.8f,%.8f,%.8f), decay=(%.8f,%.8f,%.8f)\npos=(%.8f,%.8f,%.8f), dirPdf=%.8f, oPdf=%.8f, cosine=%.8f\n" ,
			i , connectedPath[i].color.x , connectedPath[i].color.y , connectedPath[i].color.z ,
			connectedPath[i].getRadianceDecay(dist).x , connectedPath[i].getRadianceDecay(dist).y , connectedPath[i].getRadianceDecay(dist).z ,
			connectedPath[i].origin.x , connectedPath[i].origin.y , connectedPath[i].origin.z ,
			connectedPath[i].directionProb , connectedPath[i].originProb , connectedPath[i].getCosineTerm());
		*/
	}
	return vec4f(color, prob);
}

Ray BidirectionalPathTracer::link(const Path& connectedPath, int connectIndex, int i, int j)
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

float BidirectionalPathTracer::connectWeight(const Path& connectedPath, int connectIndex, vector<double>& p_forward, vector<double>& p_backward, vector<double>&distList, float expTerm)
{
	double sumExpProb = 0;

	p_forward.front() = connectedPath.front().originProb;

	for(int i=0; i<connectedPath.size()-1; i++)
	{
		p_forward[i+1] = p_forward[i] * link(connectedPath, connectIndex, i+1, i).getCosineTerm();
		distList[i] = max2((connectedPath[i+1].origin - connectedPath[i].origin).length(), EPSILON);
		
		p_forward[i+1] /= distList[i]*distList[i];
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

	p_backward.back() = connectedPath.back().originProb;

	for(int i = connectedPath.size()-1; i>0; i--)
	{
		p_backward[i-1] = p_backward[i] * link(connectedPath, connectIndex, i-1, i).getCosineTerm();
		p_backward[i-1] /= distList[i-1]*distList[i-1];
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

	for(int i=0; i<connectedPath.size()-1; i++)
	{
		if(connectedPath[i].directionSampleType == Ray::RANDOM && connectedPath[i+1].directionSampleType == Ray::RANDOM)
		{
			double p = p_forward[i] * p_backward[i+1];

			sumExpProb += pow(p, double(expTerm));
		}
	}

	if(usePT || mustUsePT(connectedPath))
	{
		sumExpProb += pow(double(p_backward.front()), double(expTerm));
	}

	double selfProb = connectIndex == -1 ? p_backward.front() : p_forward[connectIndex] * p_backward[connectIndex+1];

	if(!(pow(selfProb, double(expTerm)) / sumExpProb >= 0))
	{
		return 0;
	}

	double res = pow(selfProb, double(expTerm)) / sumExpProb;

	return res;
}

void BidirectionalPathTracer::colorByConnectingPaths(vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const Path& eyePath, const Path& lightPath, vector<unsigned>* visibilityList)
{
	unsigned maxEyePathLen = eyePath.size();
	unsigned maxLightPathLen = lightPath.size();
	unsigned maxWholePathLen = maxEyePathLen + maxLightPathLen;

	for(int wholePathLen=2; wholePathLen<=maxWholePathLen; wholePathLen++)
	{

		Path wholePath(wholePathLen);
		vector<double> p_forward(wholePathLen);
		vector<double> p_backward(wholePathLen);
		vector<double> distList(wholePathLen);

		int lplStart = (wholePathLen - int(maxEyePathLen)) > 0 ? (wholePathLen - maxEyePathLen) : 0;
		int lplEnd = wholePathLen - 1 < maxLightPathLen ? wholePathLen - 1 : maxLightPathLen;

		//if(!(wholePathLen == 3)) continue;

		for(int lightPathLen = lplStart; lightPathLen <= lplEnd; lightPathLen++)
		{
			int eyePathLen = wholePathLen - lightPathLen;

			//if (!eyePathLen == 1) continue;

			for(int li=0; li < lightPathLen; li++)
				wholePath[li] = lightPath[li];

			for(int ei=0; ei < wholePathLen - lightPathLen; ei++)
				wholePath[wholePathLen - 1 - ei] = eyePath[ei];

			int lightConnectID = lightPathLen-1;
			int eyeConnectID = lightPathLen;

			Ray& lightRay = wholePath[lightConnectID];
			Ray& eyeRay = wholePath[eyeConnectID];

			if(lightPathLen==0 && (!(wholePath.front().contactObject && wholePath.front().contactObject->emissive())))
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
		
			vec4f color_prob = connectColorProb(wholePath, lightConnectID);
			
			if(!(color_prob.w > 0))
				continue;

			if(!(vec3f(color_prob).length())>0)
				continue;

			float weight = connectWeight(wholePath, lightConnectID, p_forward, p_backward, distList);

			vec3f color = vec3f(color_prob) / color_prob.w * weight;
			
			if (eyePathLen == 1)
			{
				vec3f forward = camera.focus - camera.position;
				forward.normalize();
				vec3f dir = lightRay.origin - eyeRay.origin;
				dir.normalize();
				float cosAtCamera = abs(forward.dot(dir));
				float imageToSolidAngleFactor = camera.sightDist * camera.sightDist /
					(cosAtCamera * cosAtCamera * cosAtCamera);
				float imageToSurfaceFactor = imageToSolidAngleFactor / cosAtCamera;
				color *= (imageToSurfaceFactor / colors.size());	
				color *= cosAtCamera;
			}
			
			Ray &camRay = wholePath[wholePath.size()-1];

			if(eyePathLen > 1)
			{
				omp_set_lock(&pixelLocks[camRay.pixelID]);
				colors[camRay.pixelID] += color;
				omp_unset_lock(&pixelLocks[camRay.pixelID]);
			}
			else
			{
				vec2<float> pCoord = camera.transToPixel(wholePath[wholePath.size()-2].origin);
				int x = pCoord.x;
				int y = pCoord.y;
				if(x >= 0 && x < camera.width && y >= 0 && y < camera.height)
				{
					omp_set_lock(&pixelLocks[y*camera.width + x]);
					colors[y*camera.width + x] += color;
					omp_unset_lock(&pixelLocks[y*camera.width + x]);
				}
			}
		}
	}
}



