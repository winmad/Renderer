#include "StdAfx.h"
#include "SceneReflectiveObject.h"

Ray SceneReflectiveObject::scatter(Ray& inRay) const
{
	Ray outRay;
	vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;
	LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, position);
	vec3f normal = lf.n;
	outRay.intersectObject = NULL;
	outRay.insideObject = inRay.insideObject;
	// scatter--start
	outRay.origin = position;
	vec3f reflDir = -normal.dot(inRay.direction)*normal*2 + inRay.direction;
	reflDir.normalize();
	outRay.direction = reflDir;

	float p = maxVecComp(color);
	
	if(RandGenerator::genFloat() < p)
	{
		outRay.directionProb = p;
		outRay.contactObject = (SceneObject*) this;
		outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;
		outRay.directionSampleType = Ray::DEFINITE;

		outRay.directionSampleType = outRay.originSampleType = Ray::DEFINITE;
		outRay.color = color / outRay.getCosineTerm();
	}
	else
	{
		outRay.direction = vec3f(0, 0, 0);
		outRay.color = vec3f(0, 0, 0);
		outRay.directionProb = 1 - p;
	}
	
	// scatter--end
	return outRay;
}

float SceneReflectiveObject::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	return maxVecComp(color);
}
