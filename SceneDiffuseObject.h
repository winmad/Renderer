#pragma once
#include <algorithm>
#include "SceneObject.h"
#include "DiffuseMaterial.h"

class SceneDiffuseObject : public SceneObject
{
private:
	DiffuseMaterial *material;

public:
	DiffuseMaterial* getMaterial(){ return material; }

	SceneDiffuseObject(Scene* scene) : SceneObject(scene)
	{
		canMerge = true;
		material = new DiffuseMaterial();
		materialList.push_back(material);
	}
};

