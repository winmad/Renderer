#include "StdAfx.h"
#include "GlossyMaterial.h"

Ray GlossyMaterial::scatter(const SceneObject* object, const Ray& inRay) const
{
	Ray outRay;
	vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;
	LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, position);
	vec3f normal = lf.n;

	vec3f color = colorTex.getColor(inRay.intersectObject->getTexCoord(inRay.intersectObjectTriangleID, outRay.origin));
	vec3f coeff = coeffTex.getColor(inRay.intersectObject->getTexCoord(inRay.intersectObjectTriangleID, outRay.origin));

	outRay.insideObject = inRay.insideObject;

	// scatter--start
	outRay.origin = position;
	outRay.contactObject = (SceneObject*) object;
	outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;

	LocalFrame out_lf;
	vec3f reflDir = -normal.dot(inRay.direction)*normal*2 + inRay.direction;
	out_lf.n = reflDir;

	outRay.direction = cosineSphericalSampler.genSample(out_lf, coeff.x);

	outRay.color = bsdf.evaluate(lf, inRay.direction, outRay.direction, color, coeff.x);
	outRay.directionProb = cosineSphericalSampler.getProbDensity(out_lf, outRay.direction, coeff.x);

	if(outRay.color.length() == 0)
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

float GlossyMaterial::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	if(!outRay.contactObject)
		return 0;
	LocalFrame lf;
	vec3f normal = outRay.getContactNormal();
	vec3f reflDir = -normal.dot(inRay.direction)*normal*2 + inRay.direction;
	lf.n = reflDir;
	vec3f coeff = coeffTex.getColor(outRay.contactObject->getTexCoord(outRay.contactObjectTriangleID, outRay.origin));
	return cosineSphericalSampler.getProbDensity(lf, outRay.direction, coeff.x);
}

vec3f GlossyMaterial::getBSDF(const Ray& inRay, const Ray& outRay) const
{
	if(!outRay.contactObject)
		return vec3f(0, 0, 0);
	vec3f color = colorTex.getColor(outRay.contactObject->getTexCoord(outRay.contactObjectTriangleID, outRay.origin));
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);
	LocalFrame lf = outRay.contactObject->getAutoGenWorldLocalFrame(outRay.contactObjectTriangleID, outRay.origin);
	
	vec3f coeff = coeffTex.getColor(outRay.contactObject->getTexCoord(outRay.contactObjectTriangleID, outRay.origin));
	return p * bsdf.evaluate(lf, inRay.direction, outRay.direction, color, coeff.x);
}

float GlossyMaterial::getContinueProbability(const Ray &inRay, const Ray &outRay) const{
	vec3f color = bsdf.getColor();
	float p = *std::max_element<float*>(&color.x, (&color.x)+3);
	return p;
}

void GlossyMaterial::loadMaterialFromXML(const ConfigManager* cm, xml_node<>* node)
{
	if(node->first_node("color"))
		setColor(readVec(node->first_node("color")->value()));
	if(node->first_node("coeff"))
		setCoeff(atof(node->first_node("coeff")->value()));
	if(node->first_node("ColorTexture"))
		colorTex.loadTextureFromXML(cm, node->first_node("ColorTexture"));
	if(node->first_node("CoeffTexture"))
		coeffTex.loadTextureFromXML(cm, node->first_node("CoeffTexture"));
}