#include "StdAfx.h"
#include "AlternativeBPTracer.h"


vector<vec3f> AlternativeBPTracer::renderPixels(const Camera& camera)
{
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

		vector<unordered_map<unsigned, unsigned>> nLightRaysCaught(pixelColors.size());

		vector<Path*> lightPathList(pixelColors.size(), NULL);

#pragma omp parallel for
		for(int p=0; p<pixelColors.size(); p++)
		{
			Ray lightRay = genEmissiveSurfaceSample();
			Path *lightPath = new Path;
			samplePath(*lightPath, lightRay);
			statLightRaysCaught(nLightRaysCaught, pixelLocks, *lightPath);
			lightPathList[p] = lightPath;
		}

#pragma omp parallel for
		for(int p=0; p<pixelColors.size(); p++)
		{
			Path eyePath;

			const Path &lightPath = *lightPathList[p];

			samplePath(eyePath, camera.generateRay(p));

			if(p == pathPixelID)
				showPath = lightPath;

			colorByConnectingPaths(nLightRaysCaught, pixelLocks, renderer->camera, singleImageColors, eyePath, lightPath);

			omp_set_lock(&cmdLock);
			if(clock()-t > 100)
			{
				t = clock();
				cmd = response();
			}
			omp_unset_lock(&cmdLock);
		}

		if(cmd == "exit")
			return pixelColors;

		eliminateVignetting(singleImageColors);

		for(int i=0; i<pixelColors.size(); i++)
		{
			pixelColors[i] *= s / float(s + 1);
			pixelColors[i] += singleImageColors[i] / (s + 1)*camera.width*camera.height;
			if(clock()-t > 100)
			{
				t = clock();
				cmd = response();
			}
			delete lightPathList[i];
		}

		showCurrentResult(pixelColors);
		showPath.clear();
	}
	return pixelColors;
}

void AlternativeBPTracer::statLightRaysCaught(vector<unordered_map<unsigned, unsigned>>& nLightRaysCaught, vector<omp_lock_t> &pixelLocks, const Path& lightPath)
{
	const Camera& camera = renderer->camera;
	for(unsigned i=0; i<lightPath.size(); i++)
	{
		if(i>0 && lightPath[i].contactObject && lightPath[i].contactObject->emissive())
			break;

		if(lightPath[i].directionSampleType == Ray::DEFINITE)
			continue;

		Ray eyeRay;
		eyeRay.origin = camera.position;
		eyeRay.contactObject = (SceneObject*)&camera;
		eyeRay.insideObject = NULL;
		if(!testVisibility(link(eyeRay, lightPath[i]), link(lightPath[i], eyeRay)))
			continue;
		vec2<float> pCoord = camera.transToPixel(lightPath[i].origin);
		int x = pCoord.x;
		int y = pCoord.y;
		if(x >= 0 && x < camera.width && y >= 0 && y < camera.height)
		{
			omp_set_lock(&pixelLocks[y*camera.width + x]);
			nLightRaysCaught[y*camera.width + x][i+2] ++;
			omp_unset_lock(&pixelLocks[y*camera.width + x]);
		}
	}
}

bool AlternativeBPTracer::connectRays(Path& path, int connectIndex, int pixelID)
{
	Ray& lightRay = path[connectIndex];
	Ray& eyeRay = path[connectIndex+1];

	if(lightRay.directionSampleType == Ray::DEFINITE || eyeRay.directionSampleType == Ray::DEFINITE)
		return false;

	const Ray& prevLightRay = path[(connectIndex - 1)%path.size()];
	const Ray& prevEyeRay = path[(connectIndex + 2)%path.size()];
	lightRay.direction = eyeRay.origin - lightRay.origin;
	lightRay.direction.normalize();
	eyeRay.direction = -lightRay.direction;

	lightRay.directionProb = 1;
	lightRay.color = prevLightRay.getBSDF(lightRay);
	lightRay.directionSampleType = Ray::RANDOM;

	eyeRay.directionProb = 1;
	eyeRay.color = prevEyeRay.getBSDF(eyeRay);
	eyeRay.directionSampleType = Ray::RANDOM;

	if(connectIndex == path.size() - 2)
	{
		float connectDist = max2((lightRay.origin-eyeRay.origin).length(), EPSILON);
		float ds = connectDist*connectDist*renderer->camera.get_dw(pixelID) / lightRay.getCosineTerm();
		eyeRay.originProb = 1/ds;
	}
	else
		eyeRay.originProb = 1;
	return true;
}

vec4<float> AlternativeBPTracer::connectColorProb(const Path& connectedPath, int connectIndex)
{
	vec3f color(1, 1, 1);
	float prob = 1;

	float connectDist;

	for(unsigned i=0; i<connectedPath.size(); i++)
	{
		color *= connectedPath[i].color;

		float dist = max2((connectedPath[i+1].origin - connectedPath[i].origin).length(), EPSILON);

		if(i < connectedPath.size() - 1)
			color *= connectedPath[i].getRadianceDecay(dist);

		if(i==connectIndex && i<connectedPath.size()-1)
		{
			color *= connectedPath[i].getCosineTerm() * connectedPath[i+1].getCosineTerm() / (dist*dist);
			connectDist = dist;
		}

		if(i!=connectIndex && i!=connectIndex+1)
			color *= connectedPath[i].getCosineTerm();

		prob *= connectedPath[i].directionProb * connectedPath[i].originProb;
	}

	return vec4<float>(color, prob);
}

Ray AlternativeBPTracer::link(const Ray& src, const Ray& dst)
{
	Ray ray = src;
	ray.direction = dst.origin - src.origin;
	ray.direction.normalize();
	ray.insideObject = dst.insideObject;
	return ray;
}

float AlternativeBPTracer::connectWeight(int pixelID, int lightRaysCaught, const Path& connectedPath, int connectIndex, float expTerm)
{
	float sumExpProb = 0;

	vector<float> p_forward(connectedPath.size(), 1);
	vector<float> p_backward(connectedPath.size(), 1);
	vector<float> dist(connectedPath.size(), 0);

	p_forward.front() = connectedPath.front().originProb;

	for(int i=0; i<connectedPath.size()-1; i++)
	{
		p_forward[i+1] = p_forward[i] * link(connectedPath[i+1], connectedPath[i]).getCosineTerm();
		dist[i] = max2((connectedPath[i+1].origin - connectedPath[i].origin).length(), EPSILON);

		p_forward[i+1] /= dist[i]*dist[i];
		if(connectedPath[i].directionSampleType == Ray::RANDOM)
		{
			p_forward[i+1] *= link(connectedPath[(i-1)%connectedPath.size()], connectedPath[i]).getDirectionSampleProbDensity(link(connectedPath[i], connectedPath[i+1]));
			if(i < connectedPath.size()-2)
				p_forward[i+1] *= link(connectedPath[i], connectedPath[i+1]).getOriginSampleProbDensity(link(connectedPath[i+1], connectedPath[i+2]));
		}
	}

	p_backward.back() = connectedPath.back().originProb;

	for(int i = connectedPath.size()-1; i>0; i--)
	{
		p_backward[i-1] = p_backward[i] * link(connectedPath[i-1], connectedPath[i]).getCosineTerm() * connectedPath[i].originProb;
		p_backward[i-1] /= dist[i-1]*dist[i-1];
		if(connectedPath[i].directionSampleType == Ray::RANDOM)
		{
			p_backward[i-1] *= link(connectedPath[(i+1)%connectedPath.size()], connectedPath[i]).getDirectionSampleProbDensity(link(connectedPath[i], connectedPath[i-1]));
			if(i > 1)
				p_forward[i-1] *= link(connectedPath[i], connectedPath[i-1]).getOriginSampleProbDensity(link(connectedPath[i-1], connectedPath[i-2]));

		}
	}

	for(int i=0; i<connectedPath.size()-1; i++)
	{
		if(connectedPath[i].directionSampleType == Ray::RANDOM && connectedPath[i+1].directionSampleType == Ray::RANDOM)
		{
			float p = p_forward[i] * p_backward[i+1];
			if(i == connectedPath.size()-2)
			{
				p *= lightRaysCaught;
			}

			sumExpProb += powf(p, expTerm);
		}
	}

	sumExpProb += powf(p_backward.front(), expTerm);


	float selfProb = connectIndex == -1 ? p_backward.front() : p_forward[connectIndex] * p_backward[connectIndex+1];

	return powf(selfProb, expTerm) / sumExpProb;
}

void AlternativeBPTracer::colorByConnectingPaths(vector<unordered_map<unsigned, unsigned>>& nLightRaysCaught, vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const Path& eyePath, const Path& lightPath)
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

			if(!(eyePathLen==1)) continue;

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

			Ray &camRay = wholePath[wholePath.size()-1];
			int pixelID;
			if(eyePathLen > 1)
				pixelID = camRay.pixelID;
			else
			{
				vec2<float> pCoord = camera.transToPixel(wholePath[wholePath.size()-2].origin);
				int x = pCoord.x;
				int y = pCoord.y;
				if(x >= 0 && x < camera.width && y >= 0 && y < camera.height)
					pixelID = y*camera.width + x;
				else
					continue;
			}

			if(eyePathLen > 0 && lightPathLen > 0)
			{
				if(lightPathLen > 1 && lightRay.contactObject && lightRay.contactObject->emissive())
					continue;
				if(eyeRay.contactObject && eyeRay.contactObject->emissive())
					continue;

				if(!connectRays(wholePath, lightConnectID, pixelID))
					continue;

				if(!testVisibility(eyeRay, lightRay))
					continue;
				vec3f midPoint  = (lightRay.origin + eyeRay.origin)/2;
				eyeRay.insideObject = lightRay.insideObject = findInsideObject(midPoint, lightRay.direction);
			}

			vec4<float> color_prob = connectColorProb(wholePath, lightConnectID);

			if(color_prob.w == 0 || vec3f(color_prob).length() == 0)
				continue;
			

			int lightRays = 0;
			unordered_map<unsigned, unsigned>::iterator it = nLightRaysCaught[pixelID].find(wholePathLen);
			if(it != nLightRaysCaught[pixelID].end())
				lightRays = it->second;

			float weight = connectWeight(pixelID, lightRays, wholePath, lightConnectID);

			vec3f color = vec3f(color_prob) / color_prob.w * weight;

			omp_set_lock(&pixelLocks[pixelID]);
			colors[pixelID] += color;
			omp_unset_lock(&pixelLocks[pixelID]);
		}
	}
}
