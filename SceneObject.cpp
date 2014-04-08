#include "StdAfx.h"
#include "SceneObject.h"
#include "UniformSphericalSampler.h"

SceneObject* SceneObject::findInsideObject(const Ray& ray) const
{ 
	return scene->findInsideObject(ray, this); 
}

void SceneObject::preprocessEmissionSampler()
{
	totalArea = 0;
	for(unsigned k=0; k<getTriangleNum(); k++)
	{
		totalArea += getTriangleArea(k);
		areaValues.push_back(totalArea);
	}
	weight = totalArea;
}

Ray SceneObject::emit() const
{
	Ray ray;
	float rnd = RandGenerator::genFloat()*totalArea;
	if(!areaValues.size())
	{
		ray.direction = vec3f(0, 0, 0);
		ray.directionProb = 1;
		ray.color = vec3f(0, 0, 0);
		return ray;
	}
	unsigned index = (lower_bound(areaValues.begin(), areaValues.end(), rnd)-areaValues.begin());
	if(index == areaValues.size())
		index --; 
	ray.contactObject = (SceneObject*)this;
	ray.contactObjectTriangleID = index;
	ray.origin = genRandTrianglePosition(ray.contactObjectTriangleID);
	//UniformSphericalSampler uniformSphericalSampler;
	CosineSphericalSampler cosineSphericalSampler;

	LocalFrame lf = ray.contactObject->getAutoGenWorldLocalFrame(ray.contactObjectTriangleID, ray.origin);

	//ray.direction = uniformSphericalSampler.genSample(lf);
	ray.direction = cosineSphericalSampler.genSample(lf);

	ray.insideObject = scene->findInsideObject(ray, ray.contactObject);
	ray.current_tid = scene->getContactTreeTid(ray);
	ray.color = ray.getBSDF(ray);
	if(!emissive())
		ray.color = vec3f(1, 1, 1);

	//ray.directionProb = uniformSphericalSampler.getProbDensity(lf , ray.direction);
	ray.directionProb = cosineSphericalSampler.getProbDensity(lf, ray.direction);

	ray.originProb = weight / totalArea;
	ray.directionSampleType = ray.originSampleType = Ray::RANDOM;

	if(ray.getContactNormal().dot(ray.direction) < 0)
		ray.direction = -ray.direction;

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
	lf.n = getWorldNormal(fi, position, flat);
	vec3f axis = vec3f(0, 1, 0).cross(lf.n);
	float angle = acos(clampf(vec3f(0, 1, 0).dot(lf.n), -1, 1));
	axis.normalize();
	lf.s = vec3f(rotMat(axis, angle)*vec4<float>(vec3f(1,0,0), 0));
	lf.t = vec3f(rotMat(axis, angle)*vec4<float>(vec3f(0,0,1), 0));
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