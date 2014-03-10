#pragma once
#include "RandGenerator.h"
#include "LocalFrame.h"
#include "UniformSphericalSampler.h"
#include "macros.h"

class IsotropicPhaseSampler : public UniformSphericalSampler{
public:
	vec3f genSample(const LocalFrame &lf) const{
		vec3f genDir = UniformSphericalSampler::genSample(lf);
		return genDir;
	}

	float getProbability(const LocalFrame &lf, const vec3f &dir) const{
		return 1/(4*M_PI);
	}
};