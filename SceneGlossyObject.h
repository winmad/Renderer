#pragma once
#include <algorithm>
#include "SceneObject.h"
#include "GlossyMaterial.h"

class SceneGlossyObject : public SceneObject
{
private:
	GlossyMaterial *material;
public:
	GlossyMaterial* getMaterial(){ return material; }
	SceneGlossyObject(Scene* scene) : SceneObject(scene)
	{
		material = new GlossyMaterial();
		materialList.push_back(material);
	}
};

