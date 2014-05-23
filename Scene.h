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

	float totalArea;
	float totalVolume;
	KDTree tree;
	vector<KDTree> objKDTrees;

private:
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
		Ray genSample(bool isUniform = false) const;
		float getDirectionProbDensity(const Ray& ray) const;
		float getOriginProbDensity(const Ray& ray) const;
		void print();
	};

	struct VolumeSampler
	{
		Scene* scene;
		VolumeSampler(Scene* scene) { this->scene = scene; }
		float totalWeight;
		vector<float> weightValues;
		vector<SceneObject*> targetObjects;

		float totalVolume;

		void preprocess();
		Ray genSample(bool isUniform = false) const;
		void print();
	};

public:
	SurfaceSampler *emissiveSurfaceSampler;
	SurfaceSampler *otherSurfaceSampler;
	VolumeSampler *volumeSampler;

	Scene()
	{
		emissiveSurfaceSampler = otherSurfaceSampler = NULL;
		volumeSampler = NULL;
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
	void preprocessVolumeSampler();
	Ray genEmissionSample(bool isUniform = false) const;
	Ray genOtherSample(bool isUniform = false) const;
	Ray genVolumeSample(bool isUniform = false) const;

	float intersect(const Ray& ray, ObjSourceInformation& objSource, const KDTree::Condition* condition = NULL);
	void fillIntersectObject(vector<Ray>& rays);
	SceneObject* findInsideObject(const Ray& ray, const SceneObject* currentObject = NULL);
	bool checkInsideObject(const Ray& ray, const int insideObjectIndex);
	vector<bool> testVisibility(const vector<Ray>& rays);
	void clear();
	float getTotalArea();
	float getTotalVolume();
	float getBoundSphereRadius();
	vec3f getDiagonal();
	~Scene();
};

