#pragma once
#include "RandGenerator.h"
#include "LocalFrame.h"
#include "SphericalSampler.h"
#include "macros.h"
#include <assert.h>

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

		/** winmad's: right-hand frame */
		vec3f localFix = vec3f(sinTheta * cosPhi , cosTheta , sinTheta * sinPhi);
		vec3f res = lf.toWorld(localFix);
		return res;
		/*************/

		/** Fujun's: left-hand frame */
		vec3f localGenDir = vec3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta), u, v;
		LocalFrame newLf = generateOrthoBasis(rayDir);

		u = newLf.s;	v = newLf.t;
		
		localGenDir = u*localGenDir.x + v*localGenDir.y + rayDir*localGenDir.z;

		if ((lf.s.cross(lf.n) - lf.t).length() > 1e-6f)
			printf("error1\n");
		if ((newLf.s.cross(newLf.n) + newLf.t).length() > 1e-6f)
			printf("error2\n");
		if (abs(localGenDir.dot(rayDir) - res.dot(rayDir)) > 1e-6f)
			printf("error3\n");
		printf("========================\n");
		printf("frame1: x=(%.8f,%.8f,%.8f), y=(%.8f,%.8f,%.8f), z=(%.8f,%.8f,%.8f)\n" ,
			lf.s.x , lf.s.y , lf.s.z , lf.n.x , lf.n.y , lf.n.z , lf.t.x , lf.t.y , lf.t.z);
		printf("frame2: x=(%.8f,%.8f,%.8f), y=(%.8f,%.8f,%.8f), z=(%.8f,%.8f,%.8f)\n" ,
			newLf.s.x , newLf.s.y , newLf.s.z , newLf.n.x , newLf.n.y , newLf.n.z ,
			newLf.t.x , newLf.t.y , newLf.t.z);
		printf("(%.8f,%.8f,%.8f), (%.8f,%.8f,%.8f)\n" , res.x , res.y , res.z ,
			localGenDir.x , localGenDir.y , localGenDir.z);

		return localGenDir;
	}

	float getProbDensity(const LocalFrame& lf, const vec3f& dir) const
	{
		float cosTheta = lf.n.dot(dir);
		float temp = 1 + g*g - 2*g*cosTheta;
		
		return 1/(4*M_PI) * (1-g*g) / (temp*sqrt(temp));
	}
	
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
		lf.n = w;
		return lf;
	}
};
