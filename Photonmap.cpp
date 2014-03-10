#include "StdAfx.h"
#include "Photonmap.h"
#include "NoSelfIntersectionCondition.h"
#include <omp.h>
#include "macros.h"
//#define ENABLE_PPM

vector<vec3f> Photonmap::renderPixels(const Camera& camera)
{
	/** 
	 * prepare eye rays
	 */
	//vector<vec3f> eyeRayDirections = camera.generateRays();
	//vector<vec3f> pixelColors(eyeRayDirections.size(), vec3f(0, 0, 0));
	//vector<Ray>   eyeRays(eyeRayDirections.size());
	vector<vec3f> pixelColors(camera.width * camera.height, vec3f(0, 0, 0));
	
	preprocessEmissionSampler();
 


//#pragma omp parallel for
//	for(int i=0; i<eyeRayDirections.size(); i++)
//	{
//		Ray &eyeRay = eyeRays[i];
//		eyeRay.color = vec3f(1, 1, 1);
//		eyeRay.origin = camera.position;
//		eyeRay.direction = eyeRayDirections[i];
//		eyeRay.contactObject = NULL;
//		IntersectInfo info = intersect(eyeRay.origin, eyeRay.direction);
//		eyeRay.intersectObject = info.intersectObject;
//		eyeRay.intersectObjectTriangleID = info.triangleID;
//		eyeRay.intersectDist = info.dist;
//		eyeRay.insideObject = findInsideObject(eyeRay.origin, eyeRay.direction);
//		eyeRay.directionSampleType = Ray::DEFINITE;
//	}
 	/** 
	 * for each spp, generate an eye ray then merge photons
	 */
	spp = INT_MAX;
	std::cout << "Spp = " << spp << std::endl;


	omp_lock_t cmdLock;
	omp_init_lock(&cmdLock);
	for(int s=0; s<spp; s++)
	{
		mLightVertices.clear();
		mLightVertices.reserve(photonsWant);
		
		mVolLightVertices.clear();
		mVolLightVertices.reserve(photonsWant);
		/**
		* Shook photons for surface and volume
		*/
		shootPhotons();
		/**
		* Build global photon map & volume photon map
		*/

		mHashgrid.Reserve(mLightVertices.size());
		mHashgrid.Build(mLightVertices, mBaseRadius, photonsWant);

		mVolHashgrid.Reserve(mVolLightVertices.size());
		mVolHashgrid.Build(mVolLightVertices, mBaseRadius,photonsWant);
 

		std::cout << "sur hashgrid " << mLightVertices.size() << std::endl;
		std::cout << "vol hashgrid " << mVolLightVertices.size() << std::endl;

		vector<vec3f> singleImageColors(pixelColors.size(), vec3f(0, 0, 0));

#pragma omp parallel for
		for(int p=0; p<pixelColors.size(); p++)
		{
			pixelColors[p] *= s/float(s+1);

			Path eyePath;

			Ray eyeRaysP = camera.generateRay(p);
			eyeRaysP.directionSampleType = Ray::DEFINITE;
			eyeRaysP.originProb = eyeRaysP.directionProb = 1;
			eyeRaysP.color = vec3f(1,1,1);
			sampleMergePath(eyePath, eyeRaysP/*eyeRays[p]*/, 0);
			
			vec3f color = vec3f(0,0,0), volColor = vec3f(0,0,0);
			
			bool merge = false;
			int  mergeIndex = -1;
			vec3f eyeColor = vec3f(1,1,1);
			float eyeProb = 1.0;
			for(int i = 0; i < eyePath.size(); i++){
				if((eyePath[i].directionSampleType == Ray::RANDOM) || eyePath[i].photonType == Ray::HITVOL){
					merge = true;
					if(eyePath[i].photonType == Ray::HITVOL){
						eyeColor *= (eyePath[i].color * eyePath[i].getCosineTerm());
						eyeProb *= eyePath[i].directionProb * eyePath[i].originProb;
					}
					mergeIndex = i;
					break;
				}
				eyeColor *= eyePath[i].color;
				eyeColor *= eyePath[i].getCosineTerm();
				if(i != eyePath.size() - 1){
					float dist = max2((eyePath[i].origin - eyePath[i+1].origin).length(), EPSILON);
					eyeColor *= eyePath[i].getRadianceDecay(dist);
				}
				eyeProb *= eyePath[i].directionProb * eyePath[i].originProb;
			}
			if(!merge)	continue;
		 

			Ray keyDiffRay = eyePath[mergeIndex];
			
			bool hitVolume = keyDiffRay.photonType == Ray::HITVOL;
			if(!hitVolume){
				keyDiffRay.direction = -eyePath[mergeIndex-1].direction;
				RangeQuery query(keyDiffRay.contactObject, keyDiffRay, keyDiffRay.origin, false);

				omp_set_lock(&cmdLock);
				mHashgrid.Process(mLightVertices, query, color, false);
				omp_unset_lock(&cmdLock);

				color *= eyeColor / eyeProb;
			}
			else{
				volColor = volMerge(keyDiffRay.insideObject, keyDiffRay);
				volColor *= eyeColor / eyeProb;
			}
			//pixelColors[p] += color/(s+1) + volColor/(s+1);
			singleImageColors[p] = color/(s+1) + volColor/(s+1);
		}


		for(int i=0; i<pixelColors.size(); i++)
		{
			pixelColors[i] *= s / float(s + 1);
			pixelColors[i] += singleImageColors[i];
		}


		showCurrentResult(pixelColors);
		showPath.clear();
	}
	
	cvWaitKey(0);
	return pixelColors;

}
 
vec3f Photonmap::volMerge(SceneObject *mObj, Ray &keyDiffRay)
{
	
	omp_lock_t cmdLock;
	omp_init_lock(&cmdLock);

	
	SceneVPMObject *obj = static_cast<SceneVPMObject*>(mObj);
	float stepSize = obj->stepSize;
	float albedo = obj->getAlbedo();


	NoSelfIntersectionCondition cond(&renderer->scene, keyDiffRay);
	IntersectInfo eyeInfo = intersect(keyDiffRay.origin, keyDiffRay.direction, &cond);
	float volDist = eyeInfo.dist;
	int N = int(volDist / stepSize);
	if(N == 0){
		vec3f backGroundColor = vec3f(0,0,0);
		Ray conditionRay = keyDiffRay;
		conditionRay.origin += conditionRay.direction * volDist;
		
		NoSelfIntersectionCondition cond2(&renderer->scene, conditionRay);

		IntersectInfo backGndInfo = intersect(conditionRay.origin, conditionRay.direction, &cond2);
		if(backGndInfo.intersectObject){
			Ray backGroundMergeRay = keyDiffRay;
			backGroundMergeRay.insideObject = NULL;
			backGroundMergeRay.contactObject = backGndInfo.intersectObject;
			backGroundMergeRay.contactObjectTriangleID = backGndInfo.triangleID;
			backGroundMergeRay.origin += keyDiffRay.direction * (backGndInfo.dist);
			backGroundMergeRay.direction = -keyDiffRay.direction;
			RangeQuery query(backGndInfo.intersectObject, backGroundMergeRay, backGroundMergeRay.origin, false);
			omp_set_lock(&cmdLock);
			mHashgrid.Process(mLightVertices, query, backGroundColor, false);
			omp_unset_lock(&cmdLock);
			backGroundColor *= obj->getRadianceDecay(backGroundMergeRay, volDist);
		}
		return backGroundColor;
	}

	float step = volDist / N;
	vec3f totalVolColor(0,0,0), totalSingleColor(0,0,0), Tr(0,0,0);
	Ray outray = keyDiffRay;
	outray.direction = -keyDiffRay.direction;
	outray.contactObject = NULL;
	outray.insideObject = mObj;
	float offset = step * RandGenerator::genFloat();
	float t = offset;
	outray.origin += offset * keyDiffRay.direction;
	Tr = obj->getRadianceDecay(outray, volDist);
	for(int i = 0; i < N; i++, t+=step){
		vec3f volColor = vec3f(0,0,0);
		RangeQuery volQuery(outray.insideObject, outray, outray.origin, true);
		omp_set_lock(&cmdLock);
		mVolHashgrid.Process(mVolLightVertices, volQuery, volColor, true);
		omp_unset_lock(&cmdLock);
		
		totalVolColor += volColor * obj->getRadianceDecay(outray, t);
		//// single-scattering
		//Ray lightShadow = genEmissiveSurfaceSample();
		//vec3f pos = outray.origin;
		//vec3f connect = pos - lightShadow.origin;
		//float connectLen = connect.length();
		//connect.normalize();
		//lightShadow.direction = connect;
		//NoSelfIntersectionCondition condShadow(&renderer->scene, lightShadow);
		//IntersectInfo shadowInfo = intersect(lightShadow.origin, lightShadow.direction, &condShadow);
		//if(shadowInfo.intersectObject != obj)
		//	continue;
		//float airLen = shadowInfo.dist;
		//float volLen = connectLen - airLen;
		//vec3f lightDecay = obj->getRadianceDecay(lightShadow, volLen);
		//float lightPdf = lightShadow.getOriginSampleProbDensity(lightShadow) * lightShadow.getDirectionSampleProbDensity(lightShadow);
		//vec3f oneTimeSingleColor = obj->ds * obj->bsdf->evaluate(LocalFrame(), connect, outray.direction) * 
		//	lightShadow.color * lightDecay * obj->getRadianceDecay(outray, t) / lightPdf;
		//totalSingleColor += oneTimeSingleColor;
		
		outray.origin += keyDiffRay.direction * step;
	}

	vec3f backGroundColor = vec3f(0,0,0);
	Ray conditionRay = keyDiffRay;
	conditionRay.origin += conditionRay.direction * volDist;
	NoSelfIntersectionCondition cond2(&renderer->scene, conditionRay);
	
	IntersectInfo backGndInfo = intersect(conditionRay.origin, conditionRay.direction, &cond2);
	if(backGndInfo.intersectObject){
		Ray backGroundMergeRay = keyDiffRay;
		backGroundMergeRay.insideObject = NULL;
		backGroundMergeRay.contactObject = backGndInfo.intersectObject;
		backGroundMergeRay.contactObjectTriangleID = backGndInfo.triangleID;
		backGroundMergeRay.origin += keyDiffRay.direction * (backGndInfo.dist);
		backGroundMergeRay.direction = -keyDiffRay.direction;
		RangeQuery query(backGndInfo.intersectObject, backGroundMergeRay, backGroundMergeRay.origin, false);
		mHashgrid.Process(mLightVertices, query, backGroundColor, false);
		backGroundColor *= obj->getRadianceDecay(backGroundMergeRay, volDist);
	}
	return totalVolColor * step/*+ totalSingleColor * step*/ + backGroundColor;
}

void Photonmap::shootPhotons() 
{
	int photonsCount = 0;
	while(photonsCount < photonsWant){
 
		Path lightPath, photonPath, volPhotonPath;
		Ray lightRay = genEmissiveSurfaceSample();
		lightRay.isDirectLightPhoton = true;
		samplePath(lightPath, lightRay);
		convLightpath2Photonpath(lightPath, photonPath, volPhotonPath);
		for(int i = 0; i < photonPath.size(); i++)
			mLightVertices.push_back(photonPath[i]);
		for(int i = 0; i < volPhotonPath.size(); i++)
			mVolLightVertices.push_back(volPhotonPath[i]);
		photonsCount++;
	}
	return ;
}
 
void Photonmap::showPhotons()
{
	int size = mLightVertices.size();
	for(int i = 0; i < size; i++){
		Ray &photon = mLightVertices[i];
		std::cout << "pos " << photon.origin << " dir " << photon.direction << " f " << photon.color << " p " << photon.directionProb << std::endl;
	}
	return ;
}

void Photonmap::convLightpath2Photonpath(Path &lightPath, Path &photonPath, Path &volPhotonPath)
{
	int lightPathSize = lightPath.size();
	vec3f curC = vec3f(1,1,1);
	float curP = 1;
	
	for(int i = 0; i < lightPathSize - 1; i++){
		Ray curRay = lightPath[i], nextRay = lightPath[i+1];
		float dist = max2((curRay.origin - nextRay.origin).length(),EPSILON);
		if(i > 0 && lightPath[i].contactObject && lightPath[i].contactObject->emissive()){
			break;
		}
		
		curC *= curRay.color * curRay.getRadianceDecay(dist) * curRay.getCosineTerm();
		curP *= curRay.directionProb * curRay.originProb;
	
		Ray photon;
		photon.color = curC;
		photon.directionProb  = curP * nextRay.originProb;
		photon.direction = curRay.direction;
		photon.origin = nextRay.origin;
		photon.photonType = nextRay.photonType;
		photon.isDirectLightPhoton = curRay.isDirectLightPhoton;
		
			if(photon.photonType == Ray::OUTVOL){
				photonPath.push_back(photon);
			}
			else if(photon.photonType == Ray::INVOL){
				volPhotonPath.push_back(photon);
			}	
	}
	return ;
}

void Photonmap::sampleMergePath(Path& path, Ray& prevRay, unsigned depth) const
{
	if((prevRay.directionSampleType == Ray::RANDOM) || prevRay.photonType == Ray::HITVOL){
		path.push_back(prevRay);
		return ;
	}
	Ray termRay;
	termRay.origin = prevRay.origin;
	termRay.direction = vec3f(0, 0, 0);
	termRay.color = vec3f(0, 0, 0);
	termRay.directionProb = 1;
	termRay.insideObject = termRay.contactObject = termRay.intersectObject = NULL;
	
	termRay.directionSampleType = Ray::DEFINITE;


	path.push_back(prevRay);

	Ray nextRay;

	if(prevRay.insideObject)
	{
		nextRay = prevRay.insideObject->scatter(prevRay);
	}
	else if(prevRay.intersectObject)
	{
		nextRay = prevRay.intersectObject->scatter(prevRay);
	}
	else
	{
		path.push_back(termRay);
		return ;
	}


	if(nextRay.direction.length() < 0.5)
	{
		path.push_back(nextRay);
		return ;
	}

	if(depth+1 >= maxDepth)
	{
		path.push_back(termRay);
		return  ;
	}

	NoSelfIntersectionCondition condition(&renderer->scene, nextRay);

	Scene::ObjSourceInformation osi;
	float dist;
	dist = renderer->scene.intersect(nextRay, osi, &condition);

	if(dist < 0)
	{
		path.push_back(nextRay);
		path.push_back(termRay);
		return ;
	}
	else
	{
		nextRay.intersectObject = renderer->scene.objects[osi.objID];
		nextRay.intersectObjectTriangleID = osi.triangleID;
		nextRay.intersectDist = dist;
	}

	sampleMergePath(path, nextRay, depth + 1);
}