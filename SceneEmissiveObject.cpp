#include "StdAfx.h"
#include "SceneEmissiveObject.h"

Ray SceneEmissiveObject::scatter(const Ray& inRay , const bool russian) const
{
	Ray outRay;
	vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;
	outRay.intersectObject = NULL;
	outRay.insideObject = NULL;
	outRay.contactObject = (SceneObject*) this;
	outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;
	outRay.origin = position;
	outRay.direction = vec3f(0, 0, 0);
	outRay.color = inRay.direction.dot(outRay.getContactNormal()) <= 0 ? color : vec3f(0, 0, 0);
	outRay.directionProb = 1;
	outRay.directionSampleType = Ray::RANDOM;
	outRay.originSampleType = Ray::DEFINITE;
	return outRay;
}

vec3f SceneEmissiveObject::getBSDF(const Ray& inRay, const Ray& outRay) const
{
	if (outRay.getContactNormal().dot(outRay.direction) < 0)
		return vec3f(0.f);
	return color;
}

