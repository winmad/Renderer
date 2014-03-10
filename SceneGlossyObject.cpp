#include "StdAfx.h"
#include "SceneGlossyObject.h"
#include "CosineSphericalSampler.h"

Ray SceneGlossyObject::scatter(Ray& inRay) const
{
	Ray outRay;
	vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;
	LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, position);
	vec3f normal = lf.n;

	outRay.insideObject = inRay.insideObject;

	outRay.intersectObject = NULL;

	// scatter--start
	outRay.origin = position;
	outRay.contactObject = (SceneObject*) this;
	outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;

	CosineSphericalSampler cosineSphericalSampler;
	outRay.direction = cosineSphericalSampler.genSample(lf);

	vec3f color = bsdf.color;
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);

	if(RandGenerator::genFloat() < p)
	{
		outRay.color = bsdf.evaluate(lf, -outRay.direction, -inRay.direction);
		outRay.directionProb = p*cosineSphericalSampler.getProbDensity(lf, outRay.direction);
		outRay.photonProb = (1-p)*cosineSphericalSampler.getProbDensity(lf, outRay.direction);
	}
	else
	{
		outRay.direction = vec3f(0, 0, 0);
		outRay.color = vec3f(0, 0, 0);
		outRay.directionProb = 1;
	}

	outRay.originSampleType = Ray::DEFINITE;
	outRay.directionSampleType = Ray::RANDOM;

	// scatter--end
	return outRay;
}

float SceneGlossyObject::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	if(!outRay.contactObject)
		return 0;
	LocalFrame lf = outRay.contactObject->getAutoGenWorldLocalFrame(outRay.contactObjectTriangleID, outRay.origin);
	vec3f color = bsdf.color;
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);
	CosineSphericalSampler cosineSphericalSampler;
	return p * cosineSphericalSampler.getProbDensity(lf, outRay.direction);
}

vec3f SceneGlossyObject::getBSDF(const Ray& inRay, const Ray& outRay) const
{
	if(!outRay.contactObject)
		return vec3f(0, 0, 0);
	vec3f color = bsdf.color;
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);
	LocalFrame lf = outRay.contactObject->getAutoGenWorldLocalFrame(outRay.contactObjectTriangleID, outRay.origin);
	return p * bsdf.evaluate(lf, inRay.direction, outRay.direction);
}

float SceneGlossyObject::getContinueProbability(const Ray &inRay, const Ray &outRay) const{
	vec3f color = bsdf.color;
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);

	return p;
}
