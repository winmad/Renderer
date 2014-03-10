#pragma once
#include "DistanceSampler.h"

class HomoMediaDistSampler : public DistanceSampler
{
public:
	float sigma_t;
	HomoMediaDistSampler(float sigma_t) : sigma_t(sigma_t){}
	float sampleDist() const { 
		return -log(1.0 - RandGenerator::genFloat()) / sigma_t;
	}
	float getProbDensity(const float& dist) const { 
		return sigma_t * exp(-sigma_t * dist);
	}

};