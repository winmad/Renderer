#pragma once
#include "bsdf.h"

class GlossyBSDF : public BSDF
{
protected:
	float coeff;
	vec3f color;
public:
	void setColor(const vec3f& color) { this->color = color; }
	void setCoeff(const float& coeff) { this->coeff = coeff; }
	float getCoeff() const { return coeff; }
	vec3f getColor() const { return color; }
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay) const
	{
		if(localFrame.n.dot(inRay)<=0 && localFrame.n.dot(outRay)>=0)
		{
			vec3f reflDir = -localFrame.n.dot(inRay)*localFrame.n*2 + inRay;
			if(reflDir.dot(outRay)>0)
				return color * powf(reflDir.dot(outRay), coeff);
		}
		return vec3f(0,0,0);
	};
	virtual vec3f evaluate(const LocalFrame& localFrame, const vec3f& inRay, const vec3f& outRay, const vec3f& color, const float& coeff) const
	{
		if(localFrame.n.dot(inRay)<=0 && localFrame.n.dot(outRay)>=0)
		{
			vec3f reflDir = -localFrame.n.dot(inRay)*localFrame.n*2 + inRay;
			if(reflDir.dot(outRay)>0)
				return color * powf(reflDir.dot(outRay), coeff);
		}
		return vec3f(0,0,0);
	};
};

