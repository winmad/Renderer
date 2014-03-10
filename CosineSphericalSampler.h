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
		return RandGenerator::genHemiCosDirection(lf.n, coeff);
	}
	virtual vec3f genSample(const LocalFrame& lf, const float& coeff) const
	{
		return RandGenerator::genHemiCosDirection(lf.n, coeff);
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		return lf.n.dot(dir) > 0 ? powf(max2(lf.n.dot(dir), coeff), COS_TERM_MIN)*(coeff+1)/(2*M_PI) : 0;
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir, const float& coeff) const
	{
		return lf.n.dot(dir) > 0 ? powf(max2(lf.n.dot(dir), coeff)*(coeff+1)/(2*M_PI), COS_TERM_MIN) : 0;
	}
};