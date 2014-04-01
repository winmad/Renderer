#pragma once
#include "Ray.h"
#include "BSDF.h"
#include "smallFuncs.h"
#include "SceneObject.h"
#include "EXT/rapidxml/rapidxml.hpp"
#include "xmlHelper.h"

using namespace rapidxml;

class ConfigManager;

class Material
{
public:
	virtual float getOriginSampleProbDensity(const Ray& inRay, const Ray& outRay) { return 1; }
	virtual Ray scatter(const SceneObject* object, const Ray& inRay, const bool russian = true) const { return Ray(); }
	virtual float getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const{ return 1; }
	virtual float getContinueProbability(const Ray &inRay, const Ray &outRay) const { return 1; }
	virtual vec3f getBSDF(const Ray& inRay, const Ray& outRay) const { return vec3f(0, 0, 0); }
	virtual vec3f getRadianceDecay(const Ray& inRay, const float& dist) const { return vec3f(1, 1, 1); }
	virtual void loadMaterialFromXML(const ConfigManager* cm, xml_node<>* node){}
};

