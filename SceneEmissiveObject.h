#pragma once
#include "SceneObject.h"

class SceneEmissiveObject : public SceneObject
{
protected:
	vec3f color;
	Ray scatter(const Ray& inRay) const;
public:
	vec3f getColor() const { return color; }
	void setColor(const vec3f& color){ this->color = color; }
	SceneEmissiveObject(Scene* scene) : SceneObject(scene){
		canMerge = true;
	}
	bool emissive() const { return true; }

	vec3f getBSDF(const Ray& inRay, const Ray& outRay) const;
};

