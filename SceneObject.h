#pragma once
#include "Ray.h"
#include "simpleshape.h"
#include "Scene.h"
#include "RandGenerator.h"
#include "CosineSphericalSampler.h"
#include "NoSelfIntersectionCondition.h"
#include "Material.h"
#include "Texture.h"

class Ray;

class Scene;

class Material;

class SceneObject : public SimpleShape
{
protected:
	Scene *scene;
	vector<float> areaValues;
	vector<Material*> materialList;
public:
	float weight;
	float totalArea;
	float totalVolume;
	Texture bumpTex;
	SceneObject* findInsideObject(const Ray& ray) const;
	bool canMerge;
	virtual Ray scatter(const Ray& inRay, const bool russian = true) const;
	SceneObject(){ this->scene = NULL; }
	SceneObject(Scene* scene)
	{
		bumpTex.setColor(vec3f(1, 1, 1));
		this->scene = scene;
		unitize();
	}
	SceneObject(Scene* scene, const string &fileName, bool normalize = true) : SimpleShape(fileName, normalize)
	{
		this->scene = scene;
	}
	virtual vec3f getWorldNormal(unsigned fi, const vec3f& position, bool flat = false) const;
	virtual LocalFrame getAutoGenWorldLocalFrame(unsigned fi, const vec3f& position, bool flat = false) const;
    virtual bool hasCosineTerm(){ return true; }
	virtual bool emissive() const{ return false; }
	virtual bool isVolumeric() { return false; }
	virtual float getRefrCoeff() const{ return 1; }
	virtual void preprocessEmissionSampler();
	void normalizeEmissionWeight(float totalWeight){ weight /= totalWeight; }
	virtual Ray emit() const;
	virtual float getEmissionWeight() const { return weight; }
	virtual float getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const;
	virtual float getOriginSampleProbDensity(const Ray& inRay, const Ray& outRay) const;
	virtual float getContinueProbability(const Ray &inRay, const Ray &outRay) const;
	virtual float evaluatePdfW(const Ray &inRay, const Ray &outRay) const { return 1; }
	virtual vec3f getBSDF(const Ray& inRay, const Ray& outRay) const;
	virtual vec3f getRadianceDecay(const Ray& inRay, const float& dist) const;
	virtual ~SceneObject()
	{
		for(unsigned i=0; i<materialList.size(); i++)
			delete materialList[i];
	}
};

