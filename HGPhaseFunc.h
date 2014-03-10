#pragma once
#include "BSDF.h"

#define M_PI 3.14159265358979

class HGPhaseFunc : public BSDF
{
public:
	float g;
	HGPhaseFunc(float g) : g(g) {}
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay) const
	{
		float cosTheta = inRay.dot(outRay);
		float temp = 1 + g*g - 2*g*cosTheta;
		
		return 1/(4*M_PI) * (1-g*g) / (temp*sqrt(temp));
	};
};