#pragma once
#include "RandGenerator.h"
#include "LocalFrame.h"

class SphericalSampler
{
public:
	virtual vec3f genSample(const LocalFrame& lf) const
	{
		return vec3f(0, 0, 0);
	}
	virtual float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		return 0;
	}
};

