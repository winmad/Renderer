#pragma once
#include "sceneobject.h"
#include "RefractiveMaterial.h"

class SceneRefractiveObject : public SceneObject
{
private:
	RefractiveMaterial* material;
public:
	void setSurfColor(const vec3f& color) { material->setSurfColor(color); }
	void setDecayColor(const vec3f& color) { material->setDecayColor(color); }
	SceneRefractiveObject(Scene* scene) : SceneObject(scene)
	{ 
		material = new RefractiveMaterial;
		material->setRefrCoeff(1.5);
		material->setSurfColor(vec3f(1, 1, 1));
		material->setDecayColor(vec3f(1, 1, 1));
		materialList.push_back(material);
	}
	void setRefrCoeff(const float& coeff) { material->setRefrCoeff(coeff); }
	float getRefrCoeff() const{ return material->getRefrCoeff(); }
};


