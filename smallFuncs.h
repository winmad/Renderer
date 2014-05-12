#pragma once
#include <algorithm>
#include "macros.h"
#include "nvVector.h"
#include "nvMatrix.h"

using namespace nv;

using namespace std;

inline matrix4<float> rotMat(const vec3f& axis, const float angle)
{
	float c = cos(angle);
	float s = sin(angle);
	float _c = 1 - c;
	float _s = 1 - s;
	float x = axis.x;
	float y = axis.y;
	float z = axis.z;
	return transpose(matrix4<float>(	c+_c*x*x, _c*x*y-s*z, _c*x*z+s*y, 0,
		_c*x*y+s*z, c+_c*y*y, _c*y*z-s*x, 0,
		_c*x*z-s*y, _c*y*z+s*x, c+_c*z*z, 0,
		0, 0, 0, 1));
}

inline float maxVecComp(const vec3f& v)
{
	return *max_element((float*)&v, ((float*)&v)+3);
}

inline float clampf(float v, float min_v, float max_v)
{
	return v < min_v ? min_v : v > max_v ? max_v : v;
}

inline int clamp(int v, int min_v, int max_v)
{
	return v < min_v ? min_v : v > max_v ? max_v : v;
}

inline float y(const vec3f& v)
{
	return 0.212671f * v[0] + 0.715160f * v[1] + 0.072169f * v[2];
}

inline float intensity(const vec3f& v)
{
	return 0.212671f * v[0] + 0.715160f * v[1] + 0.072169f * v[2];
}