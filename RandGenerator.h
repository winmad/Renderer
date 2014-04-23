#pragma once
#include <time.h>
#include <random>
#include <cstdlib>
#include <stdint.h>
#include "smallFuncs.h"
#include "macros.h"

#define M_PI 3.14159265358978

using namespace nv;
using namespace std;

typedef vec3<double> vec3d;

class RNG
{
public:
	static const int N = 624;
	mutable unsigned long mt[N];
	mutable int mti;

	RNG(uint32_t _seed = 5489UL)
	{
		mti = N + 1;
		seed(_seed);
	}

	void seed(uint32_t _seed) const;
	float randFloat() const;
	uint32_t randUInt() const;
	vec3f randVector3() const;
};

class RandGenerator
{
private:
	/*
	static random_device rd;
	static tr1::mt19937 eng;
	static tr1::uniform_real_distribution<float> dis_float;
	static tr1::uniform_real_distribution<double> dis_double;
	*/
	static RNG rng;
public:
	static float genFloat();
	static vec3f genSphericalDirection();
	static vec3f genHemisphericalDirection();
	static vec3f genConcetricDisk();
	static vec3f genHemiCosDirection(float expTerm = 1, float *pdf = NULL);
};
