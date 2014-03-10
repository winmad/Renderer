#pragma once
#include "Ray.h"
#include "KDTree.h"

class Scene;

class NoSelfIntersectionCondition : public KDTree::Condition
{
public:
	Scene* scene;
	SceneObject* currentInsideObject;
	SceneObject* currentContactObject;

	NoSelfIntersectionCondition() {}
	NoSelfIntersectionCondition(Scene* scene, const Ray& ray);
	bool legal(const KDTree::Ray& ray, const KDTree::Triangle& tri, const float dist) const;
};