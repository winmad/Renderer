#pragma once
#include "LocalFrame.h"



class SceneObject;

class Ray
{
public:
	enum SampleType{DEFINITE, RANDOM};
	enum PhotonType{INVOL, OUTVOL, NOUSE, HITVOL};
	SampleType originSampleType;
	SampleType directionSampleType;
	vec3f origin;
	vec3f direction; // EndRay if length < 0.5
	vec3f color;
	float directionProb;
	float originProb;
	
	float photonProb; // TerminatingProb * oPdfW 
	PhotonType photonType;
	bool isDirectLightPhoton; 

	SceneObject* insideObject;
	SceneObject* contactObject;
	unsigned contactObjectTriangleID;

	int current_tid, intersect_tid;

	int pixelID;

// filled outside scatter()
	SceneObject* intersectObject;
	unsigned intersectObjectTriangleID;
	float intersectDist;

	Ray()
	{
		photonType = OUTVOL;
		insideObject = intersectObject = contactObject = NULL;
		color = vec3f(1, 1, 1);
		directionProb = 1;
		originProb = 1;
		photonProb = 1;
		originSampleType = DEFINITE;
		directionSampleType = RANDOM;
		pixelID = -1;
		isDirectLightPhoton = false;
		current_tid = -1;
		intersect_tid = -1;
	}

	vec3f getRadianceDecay(const float& dist) const;
	vec3f getBSDF(const Ray& outRay) const;
	float getCosineTerm(bool flat = false) const;

	float getDirectionSampleProbDensity(const Ray& outRay) const;
	float getOriginSampleProbDensity(const Ray& outRay) const;
	vec3f getContactNormal(bool flat = false) const;
};

