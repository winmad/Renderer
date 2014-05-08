#include "StdAfx.h"
#include "BidirectionalPathTracer.h"

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

			if(cmd == "exit")
				return pixelColors;

			eliminateVignetting(singleImageColors);

			for(int i=0; i<pixelColors.size(); i++)
			{
				pixelColors[i] *= s / float(s + 1);
				pixelColors[i] += singleImageColors[i] / (s + 1)*camera.width*camera.height;
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
		}

	}
	return pixelColors;
}

bool BidirectionalPathTracer::mustUsePT(const Path& connectedPath)
{
	for(unsigned i=0; i<connectedPath.size()-1; i++)
		if(connectedPath[i].directionSampleType == connectedPath[i+1].directionSampleType && connectedPath[i].directionSampleType == Ray::RANDOM)
			return false;
	return true;
}

vec4f BidirectionalPathTracer::connectColorProb(const Path& connectedPath, int connectIndex)
{
	vec3f color(1, 1, 1);
	float prob = 1;


	for(int i=0; i<connectedPath.size(); i++)
	{
		color *= connectedPath[i].color;

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

	//p_forward.front() = connectedPath.front().getOriginSampleProbDensity(connectedPath.front());
	p_forward.front() = connectedPath.front().originProb;
	//printf("%.8f , %.8f\n" , connectedPath.front().originProb , 
	//	connectedPath.front().getOriginSampleProbDensity(connectedPath.front()));

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

	//p_backward[connectedPath.size()-1] = connectedPath.back().getOriginSampleProbDensity(connectedPath.back());
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
			if(i == connectedPath.size() - 2)
				p *= renderer->camera.width * renderer->camera.height;
			sumExpProb += pow(p, double(expTerm));
		}
	}

	if(usePT || mustUsePT(connectedPath))
		sumExpProb += pow(double(p_backward.front()), double(expTerm));

	double selfProb = connectIndex == -1 ? p_backward.front() : p_forward[connectIndex] * p_backward[connectIndex+1];

	if(!(pow(selfProb, double(expTerm)) / sumExpProb >= 0))
	{
		return 0;
	}

	double res = pow(selfProb, double(expTerm)) / sumExpProb;

	// alternative
	Camera &camera = renderer->camera;
	unsigned width = camera.width, height = camera.height;

	std::vector<double> directPathProb(connectedPath.size());
	std::vector<double> reversePathProb(connectedPath.size());
	double allTechPathProb = 0;

	// connectedPath:	 Eye---------------------Light
	//           say,     V4    V3    V2    V1    V0  (direction: <-)  
	std::vector<float> distRecord;
	directPathProb.front() = connectedPath.front().originProb;//connectedPath.front().evalOriginProbability(connectedPath.front());//connectedPath.front().originProb;//
	//printf("%.8f , %.8f\n" , connectedPath.front().originProb , 
	//	connectedPath.front().evalOriginProbability(connectedPath.front()));

	for(int i = 1; i < connectedPath.size(); i++){
		float dist = MAX((connectedPath[i].origin - connectedPath[i-1].origin).length(), EPSILON);
		distRecord.push_back(dist);
		Ray linkRay = link(connectedPath, connectIndex, i, i-1);
		float cosThere = linkRay.getCosineTerm();
		directPathProb[i] = directPathProb[i-1] * cosThere / (dist * dist);
		if(connectedPath[i-1].directionSampleType == Ray::DEFINITE)
			continue;
		float linkOriProb = link(connectedPath, connectIndex, i-1, i).getOriginSampleProbDensity(
			link(connectedPath, connectIndex, i, i+1)), linkDirProb;
		if(i > 1){
			linkDirProb = link(connectedPath, connectIndex, i-2, i-1).getDirectionSampleProbDensity(
				link(connectedPath, connectIndex, i-1, i));
		}
		else{
			Ray ray = link(connectedPath, connectIndex, i-1, i);
			linkDirProb = ray.getDirectionSampleProbDensity(ray);
		}
		directPathProb[i] *= linkDirProb * linkOriProb;
	}

	// connectedPath:	 Eye---------------------Light
	//           say,     V4    V3    V2    V1    V0  (direction: ->)  
	reversePathProb.back() = connectedPath.back().originProb;//evalOriginProbability(connectedPath.back());//.originProb;
	//printf("%.8f , %.8f\n" , connectedPath.back().originProb , 
	//	connectedPath.back().evalOriginProbability(connectedPath.back()));

	for(int i = connectedPath.size()-2; i >= 0; i--){
		float dist = distRecord[i];
		Ray linkRay = link(connectedPath, connectIndex, i, i+1);
		float cosThere = linkRay.getCosineTerm();
		reversePathProb[i] = reversePathProb[i+1] * cosThere / (dist * dist);
		if(connectedPath[i+1].directionSampleType == Ray::DEFINITE)
			continue;
		float linkOriProb = link(connectedPath, connectIndex, i+1, i).getOriginSampleProbDensity(
			link(connectedPath, connectIndex, i, i-1)), linkDirProb;
		if(i < connectedPath.size()-2){
			linkDirProb = link(connectedPath, connectIndex, i+2, i-1).getDirectionSampleProbDensity(
				link(connectedPath, connectIndex, i+1, i));
		}
		else{
			Ray ray = link(connectedPath, connectIndex, i+1, i);
			linkDirProb = ray.getDirectionSampleProbDensity(ray);
		}
		reversePathProb[i] *= linkDirProb * linkOriProb;
	}

	for(int i = 0; i < connectedPath.size()-1; i++){
		if(connectedPath[i].directionSampleType==Ray::RANDOM && connectedPath[i+1].directionSampleType==Ray::RANDOM){
			double p = directPathProb[i] * reversePathProb[i+1];
			if(i == connectedPath.size()-2)
				p *= width * height;
			allTechPathProb += pow(p, double(expTerm));
		}
	}

	if(mustUsePT(connectedPath))
		allTechPathProb += pow(reversePathProb.front(), double(expTerm));

	selfProb = connectIndex == -1 ? reversePathProb.front() : directPathProb[connectIndex] * reversePathProb[connectIndex+1];

	double weight = pow(selfProb, double(expTerm)) / allTechPathProb;

	double res2 = MAX(weight, 0);
	if (abs(res - res2) > 1e-6f)
		printf("error weight\n");
	return res2;
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



