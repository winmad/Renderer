#include "StdAfx.h"
#include "RandGenerator.h"

random_device RandGenerator::rd;
tr1::mt19937 RandGenerator::eng(rd());
tr1::uniform_real_distribution<float> RandGenerator::dis_float(0, 1);
tr1::uniform_real_distribution<double> RandGenerator::dis_double(0, 1);

float RandGenerator::genFloat()
{
	return dis_float(eng);
}

vec3f RandGenerator::genSphericalDirection()
{
	float phi = acos(clampf(genFloat(), -1, 1));
	float theta = genFloat()*2*M_PI;
	vec3f dir(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
	dir = genFloat()>=0.5?dir:-dir;
	return dir;
}

vec3f RandGenerator::genHemiCosDirection(const vec3f& normal, float expTerm)
{
	float phi = acos(clampf(powf(genFloat(), 1/(expTerm+1)), -1, 1));
	float theta = genFloat()*2*M_PI;
	vec3f dir(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
	vec3f up = vec3f(0, 1, 0);
	vec3f axis = up.cross(normal);
	axis.normalize();
	float angle = acos(clampf(up.dot(normal), -1, 1));
	dir = vec3f(rotMat(axis, angle)*vec4<float>(dir, 0));
	return dir;
}

 

