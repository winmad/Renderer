#pragma once
#include "SimpleShape.h"
#include "Ray.h"
#include "SceneObject.h"
#include "KDTree.h"
#include <algorithm>

class SceneObject;

class Scene
{
public:
	struct ObjSourceInformation
	{
		unsigned objID;
		unsigned triangleID;
	};

private:
	float totalArea;
	float totalVolume;
	KDTree tree;
	vector<KDTree> objKDTrees;

	vector<unsigned> objTriangleOffsetMap;

	bool useGPU;

	struct SurfaceSampler
	{
		Scene* scene;
		SurfaceSampler(Scene* scene){ this->scene = scene; }
		float totalWeight;
		vector<float> weightValues;
		vector<SceneObject*> targetObjects;

		float totalArea;

		void preprocess();
		Ray genSample() const;
		float getDirectionProbDensity(const Ray& ray) const;
		float getOriginProbDensity(const Ray& ray) const;
	};

	SurfaceSampler *emissiveSurfaceSampler;
	SurfaceSampler *otherSurfaceSampler;
public:
	Scene()
	{
		emissiveSurfaceSampler = otherSurfaceSampler = NULL;
		//useGPU = true;
		useGPU = false;
	}
	int Scene::getContactTreeTid(const Ray& ray);
	void setGPU(bool on){ useGPU = on; }
	bool usingGPU(){ return useGPU; }
	vector<SceneObject*> objects;
	void buildKDTree();
	void buildObjKDTrees();
	void preprocessAllSamplers();
	void preprocessEmissionSampler();
	void preprocessOtherSampler();
	Ray genEmissionSample() const;
	Ray genOtherSample() const;
	float intersect(const Ray& ray, ObjSourceInformation& objSource, const KDTree::Condition* condition = NULL);
	void fillIntersectObject(vector<Ray>& rays);
	SceneObject* findInsideObject(const Ray& ray, const SceneObject* currentObject = NULL);
	vector<bool> testVisibility(const vector<Ray>& rays);
	void clear();
	float getTotalArea();
	float getTotalVolume();
	float getBoundSphereRadius();
	~Scene();
};

