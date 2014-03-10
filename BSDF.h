#pragma once
#include "nvVector.h"
#include "LocalFrame.h"

using namespace nv;

class BSDF
{
public:
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay) const
	{
		return vec3f(0, 0, 0);
	};
};

