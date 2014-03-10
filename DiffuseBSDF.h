#pragma once
#include "bsdf.h"

class DiffuseBSDF : public BSDF
{
protected:
	vec3f color;
public:
	void setColor(const vec3f& color) { this->color = color; }
	vec3f getColor() const { return color; }
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay) const
	{
		if(localFrame.n.dot(inRay)<=0 && localFrame.n.dot(outRay)>=0)
			return color/M_PI;
		else
			return vec3f(0,0,0);
	};
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay, const vec3f& color) const
	{
		if(localFrame.n.dot(inRay)<=0 && localFrame.n.dot(outRay)>=0)
			return color/M_PI;
		else
			return vec3f(0,0,0);
	};
};

