#include "StdAfx.h"
#include "Photonmap.h"
#include "NoSelfIntersectionCondition.h"
#include <omp.h>
#include "macros.h"
#include "RandGenerator.h"
#define MERGE_LEN 10
#define PPM

vector<vec3f> PhotonMap::renderPixels(const Camera& camera){
	uint width = camera.width, height = camera.height;
	std::vector<vec3f> pixelColors(width * height, vec3f(0,0,0));
		
	omp_init_lock(&surfaceHashGridLock);
	omp_init_lock(&volumeHashGridLock);
	omp_init_lock(&debugPrintLock);

	//std::vector<int> pixelMaps(pixelColors.size(), 0);

	preprocessEmissionSampler();
		
	mRadius = mBaseRadius;

	clock_t startTime = clock();

	for(uint s = 0; s < spp; s++){
		std::cout << "iteration : " << s << std::endl;
		
		std::vector<vec3f> oneIterColors(pixelColors.size(), vec3f(0,0,0));
#ifdef PPM
		if (renderer->scene.getTotalVolume() > 1e-6f)
		{
			rayMarching = true;
			mRadius = MAX(mBaseRadius * powf(powf(s + 1 , mAlpha - 1) , 1.f / 3.f) , EPSILON);
		}
		else
		{
			rayMarching = false;
			mRadius = MAX(mBaseRadius * sqrt(powf(s+1, mAlpha-1)), EPSILON);
		}
#endif
		std::vector<Path*> pixelLightPaths(pixelColors.size(), NULL);
		std::vector<LightPoint> surfaceLightVertices(0);
		std::vector<LightPoint> volumeLightVertices(0);

		surfaceHashGrid.Reserve(pixelColors.size());
		volumeHashGrid.Reserve(pixelColors.size());

#pragma omp parallel for
		// step1: sample light paths and build range search struct independently for surface and volume
		for(int p = 0; p < pixelColors.size(); p++){
			Ray lightRay = genEmissiveSurfaceSample();
			pixelLightPaths[p] = new Path;
			Path &lightPath = *pixelLightPaths[p];
			samplePath(lightPath, lightRay);
			for(int i = 1; i < lightPath.size(); i++){
				// light is not reflective
				if(lightPath[i].contactObject && lightPath[i].contactObject->emissive())
					break;
				// only store particles non-specular
				if(lightPath[i].directionSampleType == Ray::DEFINITE)
					continue;
				LightPoint lightPoint;
				lightPoint.position = lightPath[i].origin;
				lightPoint.indexInThePath = i;
				lightPoint.pathThePointIn = &lightPath;
				lightPoint.photonType = lightPath[i].photonType;
				if(lightPoint.photonType == Ray::OUTVOL){
					omp_set_lock(&surfaceHashGridLock);
					surfaceLightVertices.push_back(lightPoint);
					omp_unset_lock(&surfaceHashGridLock);
				}
				if(lightPoint.photonType == Ray::INVOL){
					omp_set_lock(&volumeHashGridLock);
					volumeLightVertices.push_back(lightPoint);
					omp_unset_lock(&volumeHashGridLock);
				}
			}
		}
		std::cout<< "vol vertices= " << volumeLightVertices.size() << " sur vertices= " << surfaceLightVertices.size() << std::endl;
			
		surfaceHashGrid.Build(surfaceLightVertices, mRadius);
		volumeHashGrid.Build(volumeLightVertices, mRadius);

		std::cout<< "finish building hashgrid" << std::endl;

		// step2: calculate pixel colors by progressive photon mapping
#pragma omp parallel for
		for(int p = 0; p < pixelColors.size(); p++){
			Path eyePath;
			if (rayMarching)
				sampleMergePath(eyePath, camera.generateRay(p), 0);
			else
				samplePath(eyePath, camera.generateRay(p));

			/*if(eyePath[1].contactObj && eyePath[1].contactObj->anisotropic()){
				pixelMaps[p] = 1;
			}*/
			throughputByDensityEstimation(oneIterColors[p], eyePath, surfaceLightVertices, volumeLightVertices);
		}
		/*std::ofstream fout(engine->renderer->name + engine->scene.name+"pixelMap.txt");
		for(int p = 0; p < pixelMaps.size(); p++)
			fout << pixelMaps[p] << ' ' ;
		fout << std::endl;
		fout.close();*/

		std::cout << "calculation done" << std::endl;

		for(uint i = 0; i < pixelColors.size(); i++){
			pixelColors[i] *= s / float(s+1);
			pixelColors[i] += camera.eliminateVignetting(oneIterColors[i], i) / (s + 1);
			delete pixelLightPaths[i];
		}
		float time = (float)(clock() - startTime) / 1000;

		if (time > recordTime)
		{
			showCurrentResult(pixelColors , &recordTime);
			recordTime += timeInterval;
		}
		showCurrentResult(pixelColors);
	}
	return pixelColors;
}

void PhotonMap::sampleMergePath(Path &path, Ray &prevRay, uint depth) const{
	path.push_back(prevRay);

	Ray terminateRay;
	terminateRay.origin = prevRay.origin;
	terminateRay.color = vec3f(0,0,0);
	terminateRay.direction = vec3f(0,0,0);
	terminateRay.directionSampleType = Ray::DEFINITE;
	terminateRay.insideObject = NULL;
	terminateRay.contactObject = NULL;
	terminateRay.intersectObject = NULL;

	Ray nextRay;
	if(prevRay.insideObject && prevRay.insideObject->isVolumeric())			
		nextRay = prevRay.insideObject->scatter(prevRay);
	else if(prevRay.intersectObject){
		if(prevRay.intersectObject->isVolumeric() && prevRay.contactObject && prevRay.contactObject->isVolumeric()){
			prevRay.origin += prevRay.direction * prevRay.intersectDist;
			prevRay.intersectDist = 0;
		}
		nextRay = prevRay.intersectObject->scatter(prevRay);
	}
	else{
		path.push_back(terminateRay);	return ;
	}

	if(nextRay.direction.length() < 0.5){
		path.push_back(nextRay);		return ;
	}
	if(depth + 1 > MERGE_LEN){
		path.push_back(terminateRay);	return ;
	}
	NoSelfIntersectionCondition condition(&renderer->scene, nextRay);
	Scene::ObjSourceInformation info;
	float dist = renderer->scene.intersect(nextRay, info, &condition);
	if(dist < 0){
		path.push_back(nextRay);
		path.push_back(terminateRay);
		return ;
	}
	else{
		nextRay.intersectObject = renderer->scene.objects[info.objID];
		nextRay.intersectObjectTriangleID = info.triangleID;
		nextRay.intersectDist = dist;
	}
	sampleMergePath(path, nextRay, depth + 1);
}

void PhotonMap::throughputByDensityEstimation(vec3f &color, Path &eyeMergePath, 
		std::vector<LightPoint> &surfaceVertices, std::vector<LightPoint> &volumeVertices)
{
	class Query{
		PhotonMap  *photonMap;
		vec3f	    contrib;
		vec3f		position;
		vec3f		hitNormal;
		float	    radius;
		int			photonsNum;
		Ray         outRay;
		float GaussianKernel(float mahalanobisDist) const{
			double exponent = exp((double)-mahalanobisDist/2);
			//photonMap->fout << " Gaussian exp = " << exponent << std::endl;
			return exponent / (2*M_PI);
		}
		float Kernel(float distSqr, float radiusSqr) const{
			float s = MAX(0, 1 - distSqr / radiusSqr);
			return 3 * s * s / M_PI;
		}
	public:
		Query(PhotonMap *map, float r, int n) : photonMap(map), radius(r), photonsNum(n)
		{}
		bool volumeMedia;
		void SetContrib(const vec3f &color) { contrib = color; }
		void SetPosition(const vec3f &pos)  { position = pos; }
		void SetOutRay(const Ray &ray) { outRay = ray; }
		void SetNormal(const vec3f &n) { hitNormal = n; }
		vec3f GetContrib() const  { return contrib; }
		vec3f GetPosition() const { return position; }
		void Process(const LightPoint &lightPoint){
			if(volumeMedia && lightPoint.photonType != Ray::INVOL)		return ;
			if(!volumeMedia && lightPoint.photonType != Ray::OUTVOL)	return ;	
			Path &lightPath = *lightPoint.pathThePointIn;
			int index = lightPoint.indexInThePath;
			vec3f photonThroughput(1,1,1);
			for(int i = 0; i < index; i++){
				photonThroughput *= lightPath[i].color / lightPath[i].directionProb / lightPath[i].originProb;
				photonThroughput *= lightPath[i].getCosineTerm();
				float dist = (lightPath[i].origin-lightPath[i+1].origin).length();
				photonThroughput *= lightPath[i].getRadianceDecay(dist);
			}
			photonThroughput /= lightPath[index].originProb;
			// runs here, photon's f/p is done.

			Ray photonRay = lightPath[index];
			photonRay.direction = lightPath[index-1].direction;
			vec3f color = photonThroughput * photonRay.getBSDF(outRay);
			float distSqr = powf((outRay.origin-lightPath[index].origin).length(), 2);
			if(intensity(color) < 1e-6f)	return ;
			float kernel = Kernel(distSqr, radius*radius);
			float normalization = volumeMedia ? kernel/(photonsNum*radius*radius*radius) : kernel/(photonsNum*radius*radius);
			//float normalization = volumeMedia==false ? 1.0 / (photonsNum*PI*radius*radius) : 1.0 / (photonsNum*PI*4.0/3*radius*radius*radius);
			contrib += color * normalization;
		}

		double sumWeight; 
		int    photonsCount;
		
		void weightScale(){
			contrib /= ( sumWeight / photonsCount );
		}

	};

	Query query(this, mRadius, mPhotonsNum);
	vec3f Tr(1,1,1), SurfaceColor(0,0,0), VolumeColor(0,0,0);
	int mergeIndex = 1;
	for(int i = 1; i < eyeMergePath.size(); i++){
		float dist = MAX((eyeMergePath[i-1].origin-eyeMergePath[i].origin).length(), EPSILON);

		if(eyeMergePath[i-1].insideObject && eyeMergePath[i-1].insideObject->isVolumeric()){
			//if(eyeMergePath[i-1].insideObject->homogeneous())
			{
				// ray marching volume radiance
				Ray volThroughRay = eyeMergePath[i-1];
				SceneVPMObject *volume = static_cast<SceneVPMObject*>(volThroughRay.insideObject);
				float stepSize = volume->stepSize;
				int N = dist / stepSize;
				if(N == 0)		N++;
				float step = dist / N;
				float offset = step * RandGenerator::genFloat();
				float t = offset;
				Tr *= volume->getRadianceDecay(volThroughRay, offset);
				for(int j = 0; j < N; j++, t+=step){
					query.SetContrib(vec3f(0,0,0));
					query.SetPosition(volThroughRay.origin + volThroughRay.direction*t);
					Ray outRay = volThroughRay;
					outRay.direction = -volThroughRay.direction;
					outRay.origin = volThroughRay.origin + volThroughRay.direction*t;
					outRay.contactObject = NULL;
					query.SetOutRay(outRay);
					query.volumeMedia = true;

					volumeHashGrid.Process(volumeVertices, query);

					Tr *= volume->getRadianceDecay(outRay, step);
					vec3f volColor = query.GetContrib();
					VolumeColor += volColor * Tr * step;
				}
			}
			/*
			else{
				// ray marching volume radiance
				
				Ray volThroughRay = eyeMergePath[i-1];
				HeterogeneousVolume *volume = static_cast<HeterogeneousVolume*>(volThroughRay.insideObj);
				float stepSize = volume->getStepSize();
				int N = dist / stepSize;
				if(N == 0)		N++;
				float step = dist / N;
				float offset = step * engine->rng->genFloat();
				float t = offset;
				Tr *= volume->randianceDecay(volThroughRay, offset);
				for(int j = 0; j < N; j++, t+=step){
					query.SetContrib(vec3f(0,0,0));
					query.SetPosition(volThroughRay.origin + volThroughRay.direction*t);
					Ray outRay = volThroughRay;
					outRay.direction = -volThroughRay.direction;
					outRay.origin = volThroughRay.origin + volThroughRay.direction*t;
					outRay.contactObj = NULL;
					query.SetOutRay(outRay);
					query.volumeMedia = true;

					volumeHashGrid.Process(volumeVertices, query);

					Tr *= volume->randianceDecay(outRay, step);
					vec3f volColor = query.GetContrib();
					VolumeColor += volColor * Tr * step;
				}
			}
			*/
		}

		if(eyeMergePath[i].contactObject && eyeMergePath[i].contactObject->emissive()){
			// eye path hit light, surface color equals to light radiance
			SurfaceColor = eyeMergePath[i].color;
			mergeIndex = i;
			break;
		}

		if(eyeMergePath[i].contactObject && eyeMergePath[i].directionSampleType == Ray::RANDOM){
			// non-specular photon density estimation
			if(eyeMergePath[i].contactObject->isVolumeric())
				continue;
			query.SetContrib(vec3f(0,0,0));
			query.SetPosition(eyeMergePath[i].origin);
			Ray outRay = eyeMergePath[i];
			outRay.direction = -eyeMergePath[i-1].direction;
			query.SetOutRay(outRay);
			query.volumeMedia = false;
			Ray fromRay = eyeMergePath[i-1];
			omp_set_lock(&surfaceHashGridLock);
			surfaceHashGrid.Process(surfaceVertices, query);
			omp_unset_lock(&surfaceHashGridLock);
			SurfaceColor = query.GetContrib();
			mergeIndex = i;
			break;
		}
	}
	color = Tr * SurfaceColor + VolumeColor;

	if (rayMarching)
	{
		for(int i = 0; i < 1/*eyeMergePath.size()-1*/; i++){
			color *= eyeMergePath[i].getCosineTerm() * eyeMergePath[i].color
				/ eyeMergePath[i].directionProb / eyeMergePath[i].originProb;
		}
	}
	else
	{
		for(int i = 0; i < mergeIndex; i++){
			color *= eyeMergePath[i].getCosineTerm() * eyeMergePath[i].color
				/ eyeMergePath[i].directionProb / eyeMergePath[i].originProb;
			if (i + 1 < mergeIndex)
			{
				float dist = (eyeMergePath[i].origin - eyeMergePath[i+1].origin).length();
				color *= eyeMergePath[i].getRadianceDecay(dist);
			}
		}
	}
}
	