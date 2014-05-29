#include "StdAfx.h"
#include "SceneObject.h"
#include "UniformSphericalSampler.h"
#include "CountHashGrid.h"

SceneObject* SceneObject::findInsideObject(const Ray& ray) const
{ 
	return scene->findInsideObject(ray, this); 
}

void SceneObject::preprocessEmissionSampler()
{
	totalArea = weight = 0.f;
	for(unsigned k=0; k<getTriangleNum(); k++)
	{
		totalArea += getTriangleArea(k);
		areaValues.push_back(totalArea);
	}
	weight = totalArea;
}

void SceneObject::preprocessOtherSampler()
{
	totalEnergy = 0.f;
	energyDensity.clear();
	energyDensity.resize(getTriangleNum());
	for (int i = 0; i < getTriangleNum(); i++)
	{
		energyDensity[i] = 0.f;
		areaValues.push_back(getTriangleArea(i));
	}
}

void SceneObject::preprocessVolumeSampler()
{
	totalVolume = volumeWeight = 0.f;
	if (!isVolumetric())
		return;
	countHashGrid = new CountHashGrid();
	countHashGrid->init(scene , objectIndex);
	countHashGrid->preprocess(scene , objectIndex);
	totalVolume = volumeWeight = countHashGrid->totVolume;
}

void SceneObject::scaleEnergyDensity(const float scale)
{
	totalEnergy *= scale;
	for (int i = 0; i < getTriangleNum(); i++)
		energyDensity[i] *= scale;
}

void SceneObject::addEnergyDensity(const int triId , const vec3f& thr)
{
	float energyDens = intensity(thr) / areaValues[triId];
	totalEnergy += energyDens;
	energyDensity[triId] += energyDens;
}

void SceneObject::singleEnergyToSumEnergy()
{
	for (int i = 1; i < getTriangleNum(); i++)
		energyDensity[i] += energyDensity[i - 1];
}

void SceneObject::sumEnergyToSingleEnergy()
{
	for (int i = getTriangleNum() - 1; i >= 1; i--)
		energyDensity[i] -= energyDensity[i - 1];
}

float SceneObject::getOriginProb(const int triId)
{
	float res;
	if (triId == 0)
		res = energyDensity[triId];
	else 
		res = energyDensity[triId] - energyDensity[triId - 1];
	res *= weight / (totalEnergy * areaValues[triId]);
	return res;
}

Ray SceneObject::emit(bool isUniformOrigin , bool isUniformDir) const
{
	Ray ray;
	
	if(!areaValues.size())
	{
		ray.direction = vec3f(0, 0, 0);
		ray.directionProb = 1;
		ray.color = vec3f(0, 0, 0);
		return ray;
	}

	if (isUniformOrigin)
	{
		float rnd = RandGenerator::genFloat()*totalArea;
		unsigned index = (lower_bound(areaValues.begin(), areaValues.end(), rnd)-areaValues.begin());
		if(index >= areaValues.size())
			index = areaValues.size() - 1; 
		ray.contactObject = (SceneObject*)this;
		ray.contactObjectTriangleID = index;
		ray.origin = genRandTrianglePosition(ray.contactObjectTriangleID);
		ray.originProb = weight / totalArea;
	}
	else
	{
		float rnd = RandGenerator::genFloat()*totalEnergy;
		unsigned index = (lower_bound(energyDensity.begin(), energyDensity.end(), rnd)-energyDensity.begin());
		if(index >= energyDensity.size())
			index = energyDensity.size() - 1; 
		ray.contactObject = (SceneObject*)this;
		ray.contactObjectTriangleID = index;
		ray.origin = genRandTrianglePosition(ray.contactObjectTriangleID);

		float prob;
		if (index == 0)
			prob = energyDensity[index] / totalEnergy;
		else
			prob = (energyDensity[index] - energyDensity[index - 1]) / totalEnergy;
		ray.originProb = weight * prob / areaValues[index];
	}

	UniformSphericalSampler uniformSphericalSampler;
	CosineSphericalSampler cosineSphericalSampler;

	LocalFrame lf = ray.contactObject->getAutoGenWorldLocalFrame(ray.contactObjectTriangleID, ray.origin);

	if (isUniformDir)
		ray.direction = uniformSphericalSampler.genSample(lf);
	else
		ray.direction = cosineSphericalSampler.genSample(lf);

	if(ray.getContactNormal().dot(ray.direction) < 0)
		ray.direction = -ray.direction;

	ray.insideObject = scene->findInsideObject(ray, ray.contactObject);

	ray.current_tid = scene->getContactTreeTid(ray);
	ray.color = ray.getBSDF(ray);
	if(!emissive())
		ray.color = vec3f(1, 1, 1);

	if (isUniformDir)
		ray.directionProb = uniformSphericalSampler.getProbDensity(lf , ray.direction) * 2.f;
	else
		ray.directionProb = cosineSphericalSampler.getProbDensity(lf, ray.direction);

	ray.directionSampleType = ray.originSampleType = Ray::RANDOM;

	if(!scene->usingGPU())
	{
		Scene::ObjSourceInformation osi;
		NoSelfIntersectionCondition condition(scene, ray);
		float dist = scene->intersect(ray, osi, &condition);
		if(dist > 0)
		{
			ray.intersectDist = dist;
			ray.intersectObject = scene->objects[osi.objID];
			ray.intersectObjectTriangleID = osi.triangleID;
		}
	}
	return ray;
}

Ray SceneObject::emitVolume(bool isUniformDir /* = false */) const
{
	Ray ray = countHashGrid->emitVolume(scene);
	//printf("%.8f , %.8f\n" , volumeWeight , ray.originProb);
	ray.originProb *= volumeWeight;
	return ray;
}

float SceneObject::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	if(&inRay == &outRay) // emitted
	{
		//UniformSphericalSampler uniformSphericalSampler;
		CosineSphericalSampler cosineSphericalSampler;

		LocalFrame lf = outRay.contactObject->getAutoGenWorldLocalFrame(outRay.contactObjectTriangleID, outRay.origin);

		//return uniformSphericalSampler.getProbDensity(lf, outRay.direction) * 2;
		return cosineSphericalSampler.getProbDensity(lf , outRay.direction);
	}
	if(materialList.size())
	{
		float prob = 0;
		for(unsigned i=0; i<materialList.size(); i++)
			prob += materialList[i]->getDirectionSampleProbDensity(inRay, outRay);
		return prob / materialList.size();
	}

	return 0;
}

vec3f SceneObject::getWorldNormal(unsigned fi, const vec3f& position, bool flat) const 
{
	vec3f original_normal = SimpleShape::getWorldNormal(fi, position, flat);
	if(bumpTex.size() <= 1 || fi >= faceVertexTexCoordIndexList.size())
		return original_normal;

	printf("use not original normal\n");

	vec3f vps[3], vts[3], vns[3];
	for(unsigned i=0; i<3; i++)
	{
		vps[i] = getWorldVertexPosition(faceVertexIndexList[fi][i]);
		if(faceVertexTexCoordIndexList[fi][i] >= vertexTexCoordList.size())
			return original_normal;
		vts[i] = vertexTexCoordList[faceVertexTexCoordIndexList[fi][i]];
	}
	vec3f uv_grad = bumpTex.getGrad(getTexCoord(fi, position));
	vec3f b1 = vps[1] - vps[0];
	vec3f b2 = vps[2] - vps[0];
	vec3f duv1 = vts[1] - vts[0];
	vec3f duv2 = vts[2] - vts[0];
	float k2 = (uv_grad.x*duv1.y - uv_grad.y*duv1.x) / (duv1.y*duv2.x - duv1.x*duv2.y);
	float k1 = (uv_grad.y - k2*duv2.y) / duv1.y;
	b1.normalize();
	b2.normalize();
	vec3f dl = k1*b1+k2*b2;
	vec3f dh = original_normal*uv_grad.z;
	if(dh.length()*1000 < dl.length())
		return original_normal;
	float angle = atan2(dh.length(), dl.length());
	vec3f axis = dl.cross(dh);
	axis.normalize();
	return vec3f(rotMat(axis, angle) * vec4f(original_normal, 0));
}

LocalFrame SceneObject::getAutoGenWorldLocalFrame(unsigned fi, const vec3f& position, bool flat) const
{
	LocalFrame lf;
	lf.buildFromNormal(getWorldNormal(fi, position, flat));
	/*
	lf.n = getWorldNormal(fi, position, flat);
	vec3f tmpT;
	tmpT = (std::abs(lf.n.z) > 0.99f) ? vec3f(1,0,0) : vec3f(0,0,1);
	lf.s = lf.n.cross(tmpT);
	lf.s.normalize();
	lf.t = lf.s.cross(lf.n);
	lf.t.normalize();
	*/
	/*
	vec3f axis = up.cross(lf.n);
	float angle = acos(clampf(up.dot(lf.n), -1, 1));
	
	if (!(axis.length() >= 1e-6f || (axis.length() < 1e-6f && n.dot(lf.n) == 1.f)))
	{
		printf("error , %.8f\n" , n.dot(lf.n));
	}
	axis.normalize();
	//if (axis.length() < 1e-6f)
	//{
	//	lf.s = vec3f(1,0,0);
	//	lf.t = vec3f(0,0,1);
	//}
	//else
	{
		lf.s = vec3f(rotMat(axis, angle)*vec4<float>(vec3f(1,0,0), 0));
		lf.t = vec3f(rotMat(axis, angle)*vec4<float>(vec3f(0,0,1), 0));
	}
	*/
	return lf;
}

float SceneObject::getOriginSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	if(&inRay == &outRay) // emitted
	{
		return weight / totalArea;
	}
	if(materialList.size())
	{
		float prob = 0;
		for(unsigned i=0; i<materialList.size(); i++)
			prob += materialList[i]->getOriginSampleProbDensity(inRay, outRay);
		return prob / materialList.size();
	}
	return 1;
}

float SceneObject::getContinueProbability(const Ray &inRay, const Ray &outRay) const
{
	return 1;
}

vec3f SceneObject::getBSDF(const Ray& inRay, const Ray& outRay) const
{
	if(materialList.size())
	{
		vec3f bsdfColor(0, 0, 0);
		for(unsigned i=0; i<materialList.size(); i++)
			bsdfColor += materialList[i]->getBSDF(inRay, outRay);
		return bsdfColor;
	}
	return vec3f(1, 1, 1);
}

vec3f SceneObject::getRadianceDecay(const Ray& inRay, const float& dist) const
{
	if(materialList.size())
	{
		vec3f decayColor(1, 1, 1);
		for(unsigned i=0; i<materialList.size(); i++)
			decayColor *= materialList[i]->getRadianceDecay(inRay, dist);
		return decayColor;
	}
	return vec3f(1, 1, 1);
}

Ray SceneObject::scatter(const Ray& inRay, const bool russian) const
{
	if(materialList.size() && this)
	{
		return materialList[rand()%materialList.size()]->scatter(this, inRay, russian);
	}
	return Ray();
}