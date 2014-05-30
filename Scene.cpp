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
}

float Scene::getTotalVolume()
{
	return volumeSampler->totalVolume;
}

void Scene::SurfaceSampler::preprocess()
{
	totalWeight = totalArea = 0.f;
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

void Scene::SurfaceSampler::preprocessForInterpath()
{
	totalWeight = totalArea = 0.f;
	weightValues.clear();
	weightValues.resize(targetObjects.size());
	for (int i = 0; i < targetObjects.size(); i++)
	{
		targetObjects[i]->preprocessOtherSampler();

		for(unsigned k=0; k<targetObjects[i]->getTriangleNum(); k++)
			totalArea += targetObjects[i]->getTriangleArea(k);
	}
}

void Scene::SurfaceSampler::normalize()
{
	totalWeight = 0.f;
	for (int i = 0; i < targetObjects.size(); i++)
	{
		totalWeight += targetObjects[i]->totalEnergy;
		weightValues[i] = totalWeight;
	}
	for (int i = 0; i < targetObjects.size(); i++)
		targetObjects[i]->weight = targetObjects[i]->totalEnergy / totalWeight;
}

void Scene::VolumeSampler::normalize()
{
	totalWeight = 0.f;
	for (int i = 0; i < targetObjects.size(); i++)
	{
		targetObjects[i]->volumeWeight = targetObjects[i]->getTotalEnergy();
		totalWeight += targetObjects[i]->volumeWeight;
		weightValues[i] = totalWeight;
	}
	for (int i = 0; i < targetObjects.size(); i++)
		targetObjects[i]->normalizeVolumeWeight(totalWeight);
}

void Scene::VolumeSampler::preprocess(bool isUniformOrigin , float mergeRadius)
{
	totalWeight = totalVolume = 0.f;
	for(unsigned i=0; i<targetObjects.size(); i++)
	{
		targetObjects[i]->preprocessVolumeSampler(isUniformOrigin , mergeRadius);
		totalWeight += targetObjects[i]->volumeWeight;
		totalVolume += targetObjects[i]->totalVolume;
		weightValues.push_back(totalWeight);
	}
	if (isUniformOrigin)
	{
		for(unsigned i=0; i<targetObjects.size(); i++)
			targetObjects[i]->normalizeVolumeWeight(totalWeight);
	}
}

Ray Scene::SurfaceSampler::genSample(bool isUniformOrigin , bool isUniformDir) const
{
	float rnd = RandGenerator::genFloat()*totalWeight;
	unsigned index = (lower_bound(weightValues.begin(), weightValues.end(), rnd)-weightValues.begin());
	if(index >= weightValues.size())
		index = weightValues.size()-1;
	return targetObjects[index]->emit(isUniformOrigin , isUniformDir);
}

Ray Scene::VolumeSampler::genSample(bool isUniformDir /* = false */) const
{
	float rnd = RandGenerator::genFloat()*totalWeight;
	unsigned index = (lower_bound(weightValues.begin(), weightValues.end(), rnd)-weightValues.begin());
	if (index >= weightValues.size())
		index = weightValues.size()-1;
	return targetObjects[index]->emitVolume(isUniformDir);
}

void Scene::SurfaceSampler::print()
{
	printf("%.8f\n" , totalWeight);
	for (int i = 0; i < targetObjects.size(); i++)
	{
		printf("%.8f %.8f\n" , targetObjects[i]->weight , targetObjects[i]->totalEnergy);
	}
}

void Scene::VolumeSampler::print()
{
	printf("%.8f\n" , totalWeight);
	for (int i = 0; i < targetObjects.size(); i++)
	{
		printf("%.8f %.8f\n" , weightValues[i] , targetObjects[i]->getTotalEnergy());
	}
}

void Scene::beginUpdateOtherSampler(const int iter)
{
	for (int i = 0; i < otherSurfaceSampler->targetObjects.size(); i++)
		otherSurfaceSampler->targetObjects[i]->sumEnergyToSingleEnergy();
	for (int i = 0; i < otherSurfaceSampler->targetObjects.size(); i++)
		otherSurfaceSampler->targetObjects[i]->scaleEnergyDensity((float)iter / ((float)iter + 1.f));
}

void Scene::beginUpdateVolumeSampler(const int iter)
{
	for (int i = 0; i < volumeSampler->targetObjects.size(); i++)
		volumeSampler->targetObjects[i]->sumEnergyToSingleEnergy();
	for (int i = 0; i < volumeSampler->targetObjects.size(); i++)
		volumeSampler->targetObjects[i]->scaleEnergyDensity((float)iter / ((float)iter + 1.f));
}

void Scene::updateOtherSampler(const int objId , const int triId , const int iter , const vec3f& thr)
{
	objects[objId]->addEnergyDensity(triId , thr / ((float)iter + 1.f));
}

void Scene::updateVolumeSampler(const int objId , const vec3f& pos , const int iter , const vec3f& thr)
{
	objects[objId]->addEnergyDensity(pos , thr / ((float)iter + 1.f));
}

void Scene::endUpdateOtherSampler()
{
	otherSurfaceSampler->normalize();
	for (int i = 0; i < otherSurfaceSampler->targetObjects.size(); i++)
		otherSurfaceSampler->targetObjects[i]->singleEnergyToSumEnergy();
}

void Scene::endUpdateVolumeSampler()
{
	volumeSampler->normalize();
	for (int i = 0; i < volumeSampler->targetObjects.size(); i++)
		volumeSampler->targetObjects[i]->singleEnergyToSumEnergy();
}

Ray Scene::genEmissionSample(bool isUniformDir) const
{
	return emissiveSurfaceSampler->genSample(true , isUniformDir);
}

Ray Scene::genOtherSample(bool isUniformOrigin , bool isUniformDir) const
{
	return otherSurfaceSampler->genSample(isUniformOrigin , isUniformDir);
}

Ray Scene::genVolumeSample(bool isUniformDir /* = false */) const
{
	return volumeSampler->genSample(isUniformDir);
}

float Scene::SurfaceSampler::getDirectionProbDensity(const Ray& ray) const
{
	//UniformSphericalSampler uniformSphericalSampler;
	CosineSphericalSampler cosineSphericalSampler;
	
	LocalFrame lf = ray.contactObject->getAutoGenWorldLocalFrame(ray.contactObjectTriangleID, ray.origin);
	
	//return uniformSphericalSampler.getProbDensity(lf, ray.direction) * 2;
	return cosineSphericalSampler.getProbDensity(lf , ray.direction);
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

void Scene::preprocessOtherSampler(bool isUniformOrigin)
{
	if(otherSurfaceSampler)
		delete otherSurfaceSampler;
	otherSurfaceSampler = new SurfaceSampler(this);
	for(unsigned i=0; i<objects.size(); i++)
	{
		if (objects[i]->canMerge && !objects[i]->isVolumetric() && !objects[i]->emissive())
			otherSurfaceSampler->targetObjects.push_back(objects[i]);
	}
	if (isUniformOrigin)
		otherSurfaceSampler->preprocess();
	else
		otherSurfaceSampler->preprocessForInterpath();
}

void Scene::preprocessVolumeSampler(bool isUniformOrigin , float mergeRadius)
{
	if(volumeSampler)
		delete volumeSampler;
	volumeSampler = new VolumeSampler(this);
	for(unsigned i=0; i<objects.size(); i++)
	{
		if (objects[i]->isVolumetric())
		{
			volumeSampler->targetObjects.push_back(objects[i]);
		}
	}
	volumeSampler->preprocess(isUniformOrigin , mergeRadius);
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
		obj->objectIndex = oi;
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

		if (obj->isVolumetric())
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
			if(currentObject == objects[i] || !objects[i]->isVolumetric())
				continue;
			dist = objKDTrees[i].intersect(kdray_back, tid);
			if(dist<1e-6f)
				continue;
			normal = objects[i]->getWorldNormal(tid, kdray_back.origin + kdray_back.direction*dist);
			if(normal.dot(kdray_back.direction) < 1e-6f)
				continue;

			dist = objKDTrees[i].intersect(kdray_front, tid);
			if(dist<1e-6f)
				continue;
			normal = objects[i]->getWorldNormal(tid, kdray_front.origin + kdray_front.direction*dist);
			if(normal.dot(kdray_front.direction) < 1e-6f)
				continue;

			return objects[i];
		}
	}

	return NULL;
}

bool Scene::checkInsideObject(const Ray& ray, const int insideObjectIndex)
{
	KDTree::Ray kdray_front, kdray_back;
	kdray_front.direction = ray.direction;
	kdray_front.origin = ray.origin;
	kdray_back.direction = -ray.direction;
	kdray_back.origin = ray.origin;

	unsigned tid;
	vec3f normal;
	float dist;
	dist = objKDTrees[insideObjectIndex].intersect(kdray_back, tid);
	if(dist<1e-6f)
		return 0;
	normal = objects[insideObjectIndex]->getWorldNormal(tid, kdray_back.origin + kdray_back.direction*dist);
	if(normal.dot(kdray_back.direction) < 1e-6f)
		return 0;

	dist = objKDTrees[insideObjectIndex].intersect(kdray_front, tid);
	if(dist<1e-6f)
		return 0;
	normal = objects[insideObjectIndex]->getWorldNormal(tid, kdray_front.origin + kdray_front.direction*dist);
	if(normal.dot(kdray_front.direction) < 1e-6f)
		return 0;

	return 1;
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