#pragma once
#include "RefractiveMaterial.h"
class SceneReflectiveObject : public SceneObject
{
private:
	RefractiveMaterial* material;
public:
	void setColor(const vec3f& color) { material->setSurfColor(color); }
	SceneReflectiveObject(Scene* scene) : SceneObject(scene)
	{
		material = new RefractiveMaterial;
		material->setDecayColor(vec3f(1, 1, 1));
		material->setSurfColor(vec3f(1, 1, 1));
		material->setRefrCoeff(-1);
		materialList.push_back(material);
		canMerge = false;
	}
};

