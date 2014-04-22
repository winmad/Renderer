// Renderer.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "SimpleShape.h"
#include "Renderer.h"

#include "IntersectionGPU.h"



int _tmain(int argc, _TCHAR* argv[])
{	
	vec3f n(0 , 0 , 1);
	CosineSphericalSampler cs;
	LocalFrame l;
	l.n = n;
	float coff = 1.f;
	vec3f d = cs.genSample(l , coff);
	fprintf(fp , "d=(%.6f , %.6f , %.6f) , pdf = %.6f\n" , d.x , d.y , d.z , cs.getProbDensity(l , d , coff));
	
	return 0;
}
