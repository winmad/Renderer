#pragma once
#include "SceneObject.h"

class RefractiveMaterial : public Material
{
protected:
	vec3f surfColor, decayColor;
	float refrCoeff;
public:
	void setRefrCoeff(const float& coeff) { refrCoeff = coeff; }
	float getRefrCoeff() const{ return refrCoeff; }
	void setSurfColor(const vec3f& color) { surfColor = color; }
	void setDecayColor(const vec3f& color) { decayColor = color; }
	virtual Ray scatter(const SceneObject* object, const Ray& inRay, const bool russian = true) const;
	virtual vec3f getRadianceDecay(const Ray& inRay, const float& dist) const { return vec3f(powf(decayColor.x, dist), powf(decayColor.y, dist), powf(decayColor.z, dist)); }
};

