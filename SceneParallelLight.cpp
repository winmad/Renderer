#include "StdAfx.h"
#include "SceneParallelLight.h"

Ray SceneParallelLight::emit()
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
	ray.contactObject = this;
	ray.contactObjectTriangleID = index;
	ray.origin = genRandTrianglePosition(ray.contactObjectTriangleID);
	LocalFrame lf = ray.contactObject->getAutoGenWorldLocalFrame(ray.contactObjectTriangleID, ray.origin);
	ray.direction = lf.n;
	if(ray.direction.dot(ray.getContactNormal()) < 0)
		ray.direction = - ray.direction;
	ray.insideObject = scene->findInsideObject(ray, ray.contactObject);
	ray.current_tid = scene->getContactTreeTid(ray);
	ray.color = ray.getBSDF(ray);
	ray.directionProb = 1;
	ray.originProb = 1 / totalArea;
	ray.directionSampleType = Ray::DEFINITE;
	ray.originSampleType = Ray::RANDOM;
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

Ray SceneParallelLight::scatter(const Ray& inRay , const bool russian) const
{
	Ray outRay;
	vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;
	outRay.intersectObject = NULL;
	outRay.insideObject = NULL;
	outRay.contactObject = (SceneObject*) this;
	outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;
	outRay.origin = position;
	outRay.direction = vec3f(0, 0, 0);
	outRay.color = vec3f(0, 0, 0);
	outRay.directionProb = 1;
	outRay.directionSampleType = Ray::DEFINITE;
	outRay.originSampleType = Ray::DEFINITE;
	return outRay;
}

float SceneParallelLight::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay)
{
	if(&inRay == &outRay) // emitted
	{
		return 1;
	}
	return 0;
}
