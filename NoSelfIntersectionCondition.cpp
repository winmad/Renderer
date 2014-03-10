#include "stdafx.h"
#include "NoSelfIntersectionCondition.h"
#include "Scene.h"

NoSelfIntersectionCondition::NoSelfIntersectionCondition(Scene* scene, const Ray& ray)
{
	this->scene = scene;
	this->currentContactObject = ray.contactObject;
	this->currentInsideObject = ray.insideObject;
}

bool NoSelfIntersectionCondition::legal(const KDTree::Ray& ray, const KDTree::Triangle& tri, const float dist) const
{
	SceneObject *intersectObject = scene->objects[((Scene::ObjSourceInformation*)tri.sourceInformation)->objID];
	unsigned fi = ((Scene::ObjSourceInformation*)tri.sourceInformation)->triangleID;
	float d = ray.direction.dot(intersectObject->getWorldNormal(fi, ray.origin + ray.direction*dist));
	bool in = d <= 0;
	bool out = d >= 0;
	if(currentInsideObject == currentContactObject)
	{
		return !(in && currentContactObject == intersectObject);
	}
	else
	{
		return !(out && currentInsideObject != intersectObject);
	}
}