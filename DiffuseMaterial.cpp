#include "StdAfx.h"
#include "DiffuseMaterial.h"

Ray DiffuseMaterial::scatter(const SceneObject* object, const Ray& inRay) const
{
	Ray outRay;
	outRay.origin = inRay.origin + inRay.direction*inRay.intersectDist;
	LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, outRay.origin);
	
	vec3f color = tex.getColor(inRay.intersectObject->getTexCoord(inRay.intersectObjectTriangleID, outRay.origin));
	
	vec3f normal = lf.n;

	outRay.insideObject = inRay.insideObject;

	outRay.intersectObject = NULL;

	// scatter--start
	outRay.contactObject = (SceneObject*) object;
	outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;

	if(inRay.direction.dot(normal)>0) // leaving
	{
		outRay.direction = inRay.direction;
		outRay.directionProb = 1;
		outRay.color = vec3f(1, 1, 1);
		outRay.photonProb = 1;
		outRay.originSampleType = Ray::DEFINITE;
		outRay.directionSampleType = Ray::DEFINITE;
		return outRay;
	}

	outRay.direction = cosineSphericalSampler.genSample(lf);

	float p = *std::max_element<float*>(&color.x, (&color.x)+3);

	if(RandGenerator::genFloat() < p)
	{
		outRay.color = bsdf.evaluate(lf, inRay.direction, outRay.direction, color);
		outRay.directionProb = p*cosineSphericalSampler.getProbDensity(lf, outRay.direction);
		outRay.photonProb = (1-p)*cosineSphericalSampler.getProbDensity(lf, outRay.direction);
	}
	else
	{
		outRay.direction = vec3f(0, 0, 0);
		outRay.color = vec3f(0, 0, 0);
		outRay.directionProb = 1 - p;
	}

	outRay.originSampleType = Ray::DEFINITE;
	outRay.directionSampleType = Ray::RANDOM;

	// scatter--end
	return outRay;
}

float DiffuseMaterial::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	if(!outRay.contactObject)
		return 0;
	LocalFrame lf = outRay.contactObject->getAutoGenWorldLocalFrame(outRay.contactObjectTriangleID, outRay.origin);
	vec3f color = tex.getColor(outRay.contactObject->getTexCoord(outRay.contactObjectTriangleID, outRay.origin));
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);
	return p * cosineSphericalSampler.getProbDensity(lf, outRay.direction);
}

float DiffuseMaterial::getContinueProbability(const Ray &inRay, const Ray &outRay) const 
{
	vec3f color = tex.getColor(outRay.contactObject->getTexCoord(outRay.contactObjectTriangleID, outRay.origin));
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);
	return p;
}

vec3f DiffuseMaterial::getBSDF(const Ray& inRay, const Ray& outRay) const 
{
	if(!outRay.contactObject)
		return vec3f(0, 0, 0);
	vec3f color = tex.getColor(outRay.contactObject->getTexCoord(outRay.contactObjectTriangleID, outRay.origin));
	LocalFrame lf = outRay.contactObject->getAutoGenWorldLocalFrame(outRay.contactObjectTriangleID, outRay.origin);
	return bsdf.evaluate(lf, inRay.direction, outRay.direction, color);
}

void DiffuseMaterial::loadMaterialFromXML(const ConfigManager* cm, xml_node<>* node)
{
	if(node->first_node("color"))
		tex.setColor(readVec(node->first_node("color")->value()));
	if(node->first_node("Texture"))
		tex.loadTextureFromXML(cm, node->first_node("Texture"));
}