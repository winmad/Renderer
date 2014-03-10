#pragma once
#include "Texture.h"
#include "SceneObject.h"
#include "GlossyBSDF.h"
#include "CosineSphericalSampler.h"

class GlossyMaterial : public Material
{
public:
protected:
	GlossyBSDF bsdf;
	CosineSphericalSampler cosineSphericalSampler;
public:
	Texture colorTex;
	Texture coeffTex;
	void setColor(const vec3f& color)
	{ 
		colorTex.setColor(color);
	}
	void setCoeff(const float& coeff)
	{
		coeffTex.setColor(coeff);
	}
	void setCoeffTex(const string& fileName) { coeffTex.loadFile(fileName); }
	void setColorTex(const string& fileName) { colorTex.loadFile(fileName); }
	virtual Ray scatter(const SceneObject* object, const Ray& inRay) const;
	virtual float getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const;
	virtual float getContinueProbability(const Ray &inRay, const Ray &outRay) const;
	virtual vec3f getBSDF(const Ray& inRay, const Ray& outRay) const;
	virtual void loadMaterialFromXML(const ConfigManager* cm, xml_node<>* node);
};

