#pragma once
#include "Texture.h"
#include "CosineSphericalSampler.h"
#include "DiffuseBSDF.h"
#include "Material.h"


class DiffuseMaterial : public Material
{
protected:
	DiffuseBSDF bsdf;
	CosineSphericalSampler cosineSphericalSampler;
public:
	Texture tex;
	virtual Ray scatter(const SceneObject* object, const Ray& inRay, const bool russian = true) const;
	virtual float getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const;
	virtual float getContinueProbability(const Ray &inRay, const Ray &outRay) const;
	virtual vec3f getBSDF(const Ray& inRay, const Ray& outRay) const;
	virtual void loadMaterialFromXML(const ConfigManager* cm, xml_node<>* node);
};

