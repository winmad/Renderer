#include "StdAfx.h"
#include "Ray.h"
#include "SceneObject.h"

vec3f Ray::getContactNormal(bool flat) const
{
	if(contactObject)
	{
		return contactObject->getWorldNormal(contactObjectTriangleID, origin, flat);
	}
	return vec3f(0, 0, 0);
}

vec3f Ray::getBSDF(const Ray& outRay) const
{
	vec3f color(0, 0, 0);
	if(outRay.contactObject)
	{
		color = outRay.contactObject->getBSDF(*this, outRay);
	}
	else if(outRay.insideObject)
	{
		color = outRay.insideObject->getBSDF(*this, outRay);
	}
	return color;
}

float Ray::getDirectionSampleProbDensity(const Ray& outRay) const
{
	float prob = 0;
	if(outRay.contactObject)
	{
		prob = outRay.contactObject->getDirectionSampleProbDensity(*this, outRay);
	}
	else if(outRay.insideObject)
	{
		prob = outRay.insideObject->getDirectionSampleProbDensity(*this, outRay);
	}
	return prob;
}

float Ray::getOriginSampleProbDensity(const Ray& outRay) const
{
	float prob = 1;

	if(insideObject)
		prob = insideObject->getOriginSampleProbDensity(*this, outRay);

	return prob;
}

float Ray::getContinueSampleProbDensity(const Ray& outRay) const
{
	float prob = 0;
	if (outRay.contactObject)
	{
		prob = outRay.contactObject->getContinueProbability(*this , outRay);
	}
	else if (outRay.insideObject)
	{
		prob = outRay.insideObject->getContinueProbability(*this , outRay);
	}
	return prob;
}

vec3f Ray::getRadianceDecay(const float& dist) const
{
	return insideObject ? insideObject->getRadianceDecay(*this, dist) : vec3f(1, 1, 1);
}

float Ray::getCosineTerm(bool flat) const
{
	float d = 1;
	if(contactObject && contactObject->hasCosineTerm() && direction.length() > 0.5)
		d = abs(getContactNormal(flat).dot(direction));
	d = max(d, COS_TERM_MIN);
	return d;
}