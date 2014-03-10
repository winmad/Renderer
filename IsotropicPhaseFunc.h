#pragma once
#include "BSDF.h"

#define M_PI 3.14159265358979

class IsotropicPhaseFunc : public BSDF
{
public:
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay) const
	{
		return 1/(4*M_PI);
	};
};