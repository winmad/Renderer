#pragma once
#include "SphericalSampler.h"

class CosineSphericalSampler : public SphericalSampler
{
private:
	float coeff;
public:
	CosineSphericalSampler(){ coeff = 1; }

	CosineSphericalSampler(const float& coeff){ this->coeff = coeff; }

	void setCoeff(const float& coeff) { this->coeff = coeff; }

	virtual vec3f genSample(const LocalFrame& lf) const
	{
		return lf.toWorld(RandGenerator::genHemiCosDirection(coeff));
	}
	virtual vec3f genSample(const LocalFrame& lf, const float& coeff) const
	{
		return lf.toWorld(RandGenerator::genHemiCosDirection(coeff));
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		float cosTheta = max2((lf.n.dot(dir)) , 0.f);
		float res = (powf(cosTheta, coeff)*(coeff+1.f)/(2.f*M_PI));
		return res;
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir, const float& coeff) const
	{
		float cosTheta = max2((lf.n.dot(dir)) , 0.f);
		float res = (powf(cosTheta, coeff)*(coeff+1.f)/(2.f*M_PI));
		return res;
	}
};