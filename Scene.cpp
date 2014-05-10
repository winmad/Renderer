#include "StdAfx.h"
#include "Scene.h"
#include "SceneEmissiveObject.h"
#include "NoSelfIntersectionCondition.h"
#include "UniformSphericalSampler.h"
#include "SceneObject.h"
#include "IntersectionGPU.h"

float Scene::getTotalArea()
{
	return otherSurfaceSampler->totalArea;
	// return 0;
	// to be filled
}

float Scene::getTotalVolume()
{
	return totalVolume;
}

void Scene::SurfaceSampler::preprocess()
{
	totalWeight = 0;
	totalArea = 0;
	for(unsigned i=0; i<targetObjects.size(); i++)
	{
		targetObjects[i]->preprocessEmissionSampler();
		totalWeight += targetObjects[i]->getEmissionWeight();
		weightValues.push_back(totalWeight);

		for(unsigned k=0; k<targetObjects[i]->getTriangleNum(); k++)
			totalArea += targetObjects[i]->getTriangleArea(k);
	}
	for(unsigned i=0; i<targetObjects.size(); i++)
		targetObjects[i]->normalizeEmissionWeight(totalWeight);
}

Ray Scene::SurfaceSampler::genSample(bool isUniform) const
{
	float rnd = RandGenerator::genFloat()*totalWeight;
	unsigned index = (lower_bound(weightValues.begin(), weightValues.end(), rnd)-weightValues.begin());
	if(index >= weightValues.size())
		index = weightValues.size()-1;
	return targetObjects[index]->emit(isUniform);
}

Ray Scene::genEmissionSample(bool isUniform) const
{
	return emissiveSurfaceSampler->genSample();
}

Ray Scene::genOtherSample(bool isUniform) const
{
	return otherSurfaceSampler->genSample();
}

float Scene::SurfaceSampler::getDirectionProbDensity(const Ray& ray) const
{
	//UniformSphericalSampler uniformSphericalSampler;
	CosineSphericalSampler cosineSphericalSampler;
	
	LocalFrame lf = ray.contactObject->getAutoGenWorldLocalFrame(ray.contactObjectTriangleID, ray.origin);
	
	//return uniformSphericalSampler.getProbDensity(lf, ray.direction) * 2;
	return cosineSphericalSampler.getProbDensity(lf , ray.direction);
}

void Scene::preprocessAllSamplers()
{
	if(emissiveSurfaceSampler)
		delete emissiveSurfaceSampler;
	if(otherSurfaceSampler)
		delete otherSurfaceSampler;
	emissiveSurfaceSampler = new SurfaceSampler(this);
	otherSurfaceSampler = new SurfaceSampler(this);
	for(unsigned i=0; i<objects.size(); i++)
	{
		if(objects[i]->emissive())
			emissiveSurfaceSampler->targetObjects.push_back(objects[i]);
		else
		{
			if (objects[i]->canMerge)
				otherSurfaceSampler->targetObjects.push_back(objects[i]);
		}
	}
	emissiveSurfaceSampler->preprocess();
	otherSurfaceSampler->preprocess();
}

void Scene::preprocessEmissionSampler()
{
	if(emissiveSurfaceSampler)
		delete emissiveSurfaceSampler;
	emissiveSurfaceSampler = new SurfaceSampler(this);
	for(unsigned i=0; i<objects.size(); i++)
	{
		if(objects[i]->emissive())
			emissiveSurfaceSampler->targetObjects.push_back(objects[i]);
	}
	emissiveSurfaceSampler->preprocess();
}

void Scene::preprocessOtherSampler()
{
	if(otherSurfaceSampler)
		delete otherSurfaceSampler;
	otherSurfaceSampler = new SurfaceSampler(this);
	for(unsigned i=0; i<objects.size(); i++)
	{
		if (objects[i]->canMerge)
			otherSurfaceSampler->targetObjects.push_back(objects[i]);
	}
	otherSurfaceSampler->preprocess();
}

vec3f Scene::getDiagonal()
{
	return tree.getDiagonal();
}

float Scene::getBoundSphereRadius()
{
	vec3f diameter = tree.getDiagonal();
	return diameter.length() * 0.5f;
}

void Scene::buildKDTree()
{
	unsigned vi_offset = 0;
	unsigned tri_offset = 0;
	for(unsigned oi=0; oi<objects.size(); oi++)
	{
		SceneObject* obj = objects[oi];
		for(unsigned i=0; i<obj->getVertexNum(); i++)
		{
			tree.vertexPositionList.push_back(obj->getWorldVertexPosition(i));
		}

		for(unsigned i=0; i<obj->getTriangleNum(); i++)
		{
			KDTree::Triangle tri;
			tri.vertexIndices = obj->getVertexIndices(i) + vec3ui(vi_offset, vi_offset, vi_offset);
			tri.sourceInformation = new ObjSourceInformation;
			((ObjSourceInformation*)tri.sourceInformation)->objID = oi;
			((ObjSourceInformation*)tri.sourceInformation)->triangleID = i;
			tree.triangleList.push_back(tri);
		}
		objTriangleOffsetMap.push_back(tri_offset);
		tri_offset += obj->getTriangleNum();
		vi_offset += obj->getVertexNum();
	}
	
	tree.build();
	if(useGPU)
	{
		glutExit();

		IntersectionGPU::init();
		
		IntersectionGPU::loadKDTreeToShader(tree);
		
	}
}

void Scene::buildObjKDTrees()
{
	objKDTrees.resize(objects.size());
	for(unsigned oi=0; oi<objects.size(); oi++)
	{
		SceneObject* obj = objects[oi];
		for(unsigned i=0; i<obj->getVertexNum(); i++)
		{
			objKDTrees[oi].vertexPositionList.push_back(obj->getWorldVertexPosition(i));
		}

		for(unsigned i=0; i<obj->getTriangleNum(); i++)
		{
			KDTree::Triangle tri;
			tri.vertexIndices = obj->getVertexIndices(i);
			tri.sourceInformation = new ObjSourceInformation;
			((ObjSourceInformation*)tri.sourceInformation)->objID = oi;
			((ObjSourceInformation*)tri.sourceInformation)->triangleID = i;
			objKDTrees[oi].triangleList.push_back(tri);
		}
		objKDTrees[oi].build();

		if (obj->isVolumeric())
		{
			float r = objKDTrees[oi].getDiagonal().length() / 2;
			obj->totalVolume = 4.f / 3.f * r * r * r;
		}
		else
		{
			obj->totalVolume = 0.f;
		}
		totalVolume += obj->totalVolume;
	}
}

void Scene::clear()
{
	for(unsigned i=0; i<objects.size(); i++)
		delete objects[i];
	objects.clear();
	for(unsigned k=0; k<objKDTrees.size(); k++)
		objKDTrees[k].destroy();
	tree.destroy();
	if(emissiveSurfaceSampler)
		delete emissiveSurfaceSampler;
	if(otherSurfaceSampler)
		delete otherSurfaceSampler;
}

Scene::~Scene()
{
	clear();
}

float Scene::intersect(const Ray& ray, ObjSourceInformation& objSource, const KDTree::Condition* condition)
{
	KDTree::Ray kdray;
	kdray.direction = ray.direction;
	kdray.origin = ray.origin;

	unsigned tid;
	float dist = tree.intersect(kdray, tid, condition);
	if(dist >= 0)
		objSource = *(ObjSourceInformation*)(tree.triangleList[tid].sourceInformation);
	return dist;
}

SceneObject* Scene::findInsideObject(const Ray& ray, const SceneObject* currentObject)
{
	if(!useGPU || true)
	{
		KDTree::Ray kdray_front, kdray_back;
		kdray_front.direction = ray.direction;
		kdray_front.origin = ray.origin;
		kdray_back.direction = -ray.direction;
		kdray_back.origin = ray.origin;

		for(int i=objKDTrees.size()-1; i>=0; i--)
		{
			unsigned tid;
			vec3f normal;
			float dist;
			if(currentObject == objects[i])
				continue;
			dist = objKDTrees[i].intersect(kdray_back, tid);
			if(dist<0)
				continue;
			normal = objects[i]->getWorldNormal(tid, kdray_back.origin + kdray_back.direction*dist);
			if(normal.dot(kdray_back.direction) < 0)
				continue;

			dist = objKDTrees[i].intersect(kdray_front, tid);
			if(dist<0)
				continue;
			normal = objects[i]->getWorldNormal(tid, kdray_front.origin + kdray_front.direction*dist);
			if(normal.dot(kdray_front.direction) < 0)
				continue;

			return objects[i];
		}
	}

	return NULL;
}

int Scene::getContactTreeTid(const Ray& ray)
{
	if(ray.contactObject)
	{
		for(int i=0; i<objects.size(); i++)
			if(ray.contactObject == objects[i])
				return objTriangleOffsetMap[i] + ray.contactObjectTriangleID;
	}
	return -1;
}

void Scene::fillIntersectObject(vector<Ray>& rays)
{
	vector<IntersectionGPU::Ray> raysInGPU(rays.size());
	for(unsigned i=0; i<rays.size(); i++)
	{
		raysInGPU[i].direction = rays[i].direction;
		raysInGPU[i].origin = rays[i].origin;
		raysInGPU[i].last_tid = rays[i].contactObject ? rays[i].current_tid : -1;
	}
	vector<vec4f> dist_tid = IntersectionGPU::query(raysInGPU);
	for(unsigned i=0; i<rays.size(); i++)
	{
		if(dist_tid[i].x >= 0)
		{
			int tid = int(dist_tid[i].y);
			ObjSourceInformation *osi = (ObjSourceInformation *)tree.triangleList[tid].sourceInformation;
			rays[i].intersectObject = objects[osi->objID];
			rays[i].intersectDist = dist_tid[i].x;
			rays[i].intersectObjectTriangleID = osi->triangleID;
			rays[i].intersect_tid = tid;
		}
	}
}

vector<bool> Scene::testVisibility(const vector<Ray>& rays)
{
	vector<bool> visibilityList(rays.size());
	vector<IntersectionGPU::Ray> raysInGPU(rays.size());
	for(unsigned i=0; i<rays.size(); i++)
	{
		raysInGPU[i].direction = rays[i].direction;
		raysInGPU[i].origin = rays[i].origin;
		raysInGPU[i].last_tid = rays[i].current_tid;
	}

	vector<vec4f> dist_tid = IntersectionGPU::query(raysInGPU);
	vector<float> shiftDists;
	vector<unsigned> visIDs;
	raysInGPU.clear();
	for(unsigned i=0; i<rays.size(); i++)
	{
		if(dist_tid[i].x < 0 || dist_tid[i].x - rays[i].intersectDist >= -EPSILON)
			visibilityList[i] = true;
		else
		{
			int tid = int(dist_tid[i].y);
			ObjSourceInformation *osi = (ObjSourceInformation *)tree.triangleList[tid].sourceInformation;
			bool legal = objects[osi->objID] != rays[i].contactObject;
			vec3f normal = objects[osi->objID]->SimpleShape::getWorldNormal(osi->triangleID, rays[i].origin+dist_tid[i].x*rays[i].direction);
			bool in = normal.dot(rays[i].direction) < 0;
			legal |= objects[osi->objID] == rays[i].insideObject && !in;
			visibilityList[i] = false;
			if(!legal)
			{
				IntersectionGPU::Ray retestRay;
				retestRay.direction = rays[i].direction;
				retestRay.origin = rays[i].origin + rays[i].direction * dist_tid[i].x;
				retestRay.last_tid = tid;
				visIDs.push_back(i);
				raysInGPU.push_back(retestRay);
				shiftDists.push_back(dist_tid[i].x);
			}
		}
	}

	// retest illegal rays
	dist_tid = IntersectionGPU::query(raysInGPU);

	for(unsigned i=0; i<raysInGPU.size(); i++)
	{
		if(dist_tid[i].x < 0 || dist_tid[i].x - rays[visIDs[i]].intersectDist + shiftDists[i] >= -EPSILON)
		{
			visibilityList[visIDs[i]] = true;
		}
	}

	return visibilityList;
}