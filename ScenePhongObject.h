#pragma once
#include <algorithm>
#include "SceneObject.h"
#include "DiffuseMaterial.h"
#include "GlossyMaterial.h"

class ScenePhongObject : public SceneObject
{
private:
	GlossyMaterial *glossyMaterial;
	DiffuseMaterial *diffuseMaterial;
public:
	GlossyMaterial* getGlossyMaterial(){ return glossyMaterial; }
	DiffuseMaterial* getDiffuseMaterial(){ return diffuseMaterial; }
	ScenePhongObject(Scene* scene) : SceneObject(scene)
	{
		glossyMaterial = new GlossyMaterial();
		materialList.push_back(glossyMaterial);
		diffuseMaterial = new DiffuseMaterial();
		materialList.push_back(diffuseMaterial);
	}
};