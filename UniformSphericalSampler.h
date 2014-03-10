#pragma once
#include "SphericalSampler.h"

class UniformSphericalSampler : public SphericalSampler
{
public:
	virtual vec3f genSample(const LocalFrame& lf) const
	{
		return RandGenerator::genSphericalDirection();
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		return 1/(4*M_PI);
	}
};