#include "stdafx.h"
#include "xmlHelper.h"

matrix4<float> readMatrix(const char* value)
{
	float v[16];

	sscanf_s(value, "[%f, %f, %f, %f][%f, %f, %f, %f][%f, %f, %f, %f][%f, %f, %f, %f]",
		v, v+1, v+2, v+3,
		v+4,v+5, v+6, v+7,
		v+8, v+9, v+10, v+11,
		v+12, v+13, v+14, v+15);
	matrix4<float> mat;
	mat.set_value(v);
	return transpose(mat);
}

vec3f readVec(const char* value)
{
	vec3f vec;
	sscanf_s(value, "(%f, %f, %f)", &vec.x, &vec.y, &vec.z);
	return vec;
}