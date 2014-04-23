#pragma once
#include "nvVector.h"

using namespace nv;



struct LocalFrame
{
	vec3f s;
	vec3f t;
	vec3f n;

	void buildFromNormal(const vec3f& normal)
	{
		n = normal;
		vec3f tmpT;
		tmpT = (std::abs(n.z) > 0.99f) ? vec3f(1,0,0) : vec3f(0,0,1);
		s = n.cross(tmpT);
		s.normalize();
		t = s.cross(n);
		t.normalize();
	}

	vec3f toWorld(const vec3f& a) const
	{
		return s * a.x + n * a.y + t * a.z;
	}

	vec3f toLocal(const vec3f &a) const
	{
		return vec3f(a.dot(s) , a.dot(n) , a.dot(t));
	}
};



//	   (n)
//		A	
//		|
//	  /-|----/
//	 /	|---->(s)
//	/--/---/
//	  /	
//	 v
//	(t)	