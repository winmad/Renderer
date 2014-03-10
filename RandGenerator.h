#pragma once
#include <time.h>
#include <random>
#include "smallFuncs.h"
#include "macros.h"

#define M_PI 3.14159265358978

using namespace nv;
using namespace std;

typedef vec3<double> vec3d;

class RandGenerator
{
private:
	static random_device rd;
	static tr1::mt19937 eng;
	static tr1::uniform_real_distribution<float> dis_float;
	static tr1::uniform_real_distribution<double> dis_double;
public:
	static float genFloat();
	static vec3f genSphericalDirection();
	static vec3f genHemiCosDirection(const vec3f& normal, float expTerm = 1);
 
};

