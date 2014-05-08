#pragma once
#include "DistanceSampler.h"

class HomoMediaDistSampler : public DistanceSampler
{
public:
	vec3f sigma_t;
	HomoMediaDistSampler(vec3f _sigma_t) : sigma_t(_sigma_t){}
	float sampleDist() const { 
		return -log(1.0 - RandGenerator::genFloat()) / y(sigma_t);
	}
	float getProbDensity(const float& dist) const { 
		return y(sigma_t) * exp(-y(sigma_t) * dist);
	}

};