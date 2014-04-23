#pragma once
#include "RandGenerator.h"
#include "LocalFrame.h"
#include "SphericalSampler.h"
#include "macros.h"

class HGPhaseSampler : public SphericalSampler{
public:
	float g;
	HGPhaseSampler(float g) : g(g) {}

	vec3f genSample(const LocalFrame& lf) const
	{
		vec3f rayDir = lf.n;	
		rayDir.normalize();

		float e1 = RandGenerator::genFloat(), e2 = RandGenerator::genFloat();
		float cosTheta;
		if(abs(g) < EPSILON){
			cosTheta = 1 - 2*e1;
		}
		else{
			float sqrTerm = (1 - g*g) / (1 - g + 2*g*e1);
			cosTheta = (1 + g*g - sqrTerm * sqrTerm) / (2*g);
		}
		float sinTheta = sqrt(1 - cosTheta * cosTheta), sinPhi, cosPhi;
		sinPhi = sin(2.0 * M_PI * e2);
		cosPhi = cos(2.0 * M_PI * e2);
		
		vec3f localGenDir = vec3f(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi), u, v;
		//LocalFrame newLf = generateOrthoBasis(rayDir);
		LocalFrame newLf;
		newLf.buildFromNormal(rayDir);
		return newLf.toWorld(localGenDir);
		/*
		u = newLf.s;	v = newLf.t;
		
		localGenDir = u*localGenDir.x + v*localGenDir.y + rayDir*localGenDir.z;
		return localGenDir;
		*/
	}

	float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		float cosTheta = lf.n.dot(dir);
		float temp = 1 + g*g - 2*g*cosTheta;
		
		return 1/(4*M_PI) * (1-g*g) / (temp*sqrt(temp));
	}
	/*
	LocalFrame generateOrthoBasis(vec3f w) const{
		vec3f coVec = w, u, v;
		if (fabs(w.x) <= fabs(w.y))
			if (fabs(w.x) <= fabs(w.z)) coVec = vec3f(0,-w.z,w.y);
			else coVec = vec3f(-w.y,w.x,0);
		else if (fabs(w.y) <= fabs(w.z)) coVec = vec3f(-w.z,0,w.x);
		else coVec = vec3f(-w.y,w.x,0);
		coVec.normalize();
		u = w.cross(coVec);
		v = w.cross(u);
		LocalFrame lf;
		lf.s = u;
		lf.t = v;
		return lf;
	}
	*/
};
