#pragma once
#include "RandGenerator.h"

class DistanceSampler
{
public:
	virtual float sampleDist() const { return 0; }
	virtual float getProbDensity(const float& dist) const { return 1; }
};