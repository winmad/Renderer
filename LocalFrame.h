#pragma once
#include "nvVector.h"

using namespace nv;



struct LocalFrame
{
	vec3f s;
	vec3f t;
	vec3f n;
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