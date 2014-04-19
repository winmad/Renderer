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
		float pdf;
		return RandGenerator::genHemiCosDirection(lf.n, coeff, &pdf);
	}
	virtual vec3f genSample(const LocalFrame& lf, const float& coeff) const
	{
		float pdf;
		return RandGenerator::genHemiCosDirection(lf.n, coeff, &pdf);
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		float res = (lf.n.dot(dir) > 1e-6f ? (powf(lf.n.dot(dir), coeff)*(coeff+1.f)/(2.f*M_PI)) : 0);
		if (lf.n.dot(dir) > 1e-6f && abs(res - lf.n.dot(dir) / M_PI) > 1e-6f)
		{
			printf("error! %.6f,%6f,%.6f\n" , coeff , res , lf.n.dot(dir) / M_PI);
		}
		return max2(res , 0.f);
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir, const float& coeff) const
	{
		float res = (lf.n.dot(dir) > 1e-6f ? (powf(lf.n.dot(dir), coeff)*(coeff+1.f)/(2.f*M_PI)) : 0);
		if (lf.n.dot(dir) > 1e-6f && abs(res - lf.n.dot(dir) / M_PI) > 1e-6f)
		{
			printf("error! %.6f,%6f,%.6f\n" , coeff , res , lf.n.dot(dir) / M_PI);
		}
		return max2(res , 0.f);
	}
};