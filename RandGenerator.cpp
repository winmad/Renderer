#include "StdAfx.h"
#include "RandGenerator.h"

#define M 397
#define MATRIX_A 0x9908b0dfUL
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

RNG RandGenerator::rng(5489UL);

void RNG::seed(uint32_t _seed) const
{
	mt[0] = _seed & 0xffffffffUL;
	for (mti = 1; mti < N; mti++)
	{
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		mt[mti] &= 0xffffffffUL;
	}
}

float RNG::randFloat() const
{
	float v = (randUInt() & 0xffffff) / float(1 << 24);
	return v;
}

uint32_t RNG::randUInt() const
{
	unsigned long y;
	static unsigned long mag01[2] = {0x0UL , MATRIX_A};

	if (mti >= N)
	{
		int k;
		if (mti == N + 1)
		{
			seed(5489UL);
		}

		for (k = 0; k < N - M; k++)
		{
			y = (mt[k] & UPPER_MASK) | (mt[k + 1] & LOWER_MASK);
			mt[k] = mt[k + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}

		for (; k < N - 1; k++)
		{
			y = (mt[k] & UPPER_MASK) | (mt[k + 1] & LOWER_MASK);
			mt[k] = mt[k + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}

		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}
	y = mt[mti++];

	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

vec3f RNG::randVector3() const
{
	float a = randFloat();
	float b = randFloat();
	float c = randFloat();
	return vec3f(a , b , c);
}

float RandGenerator::genFloat()
{
	return rng.randFloat();
}

vec3f RandGenerator::genSphericalDirection()
{
	/*
	float phi = acos(clampf(genFloat(), -1, 1));
	float theta = genFloat()*2*M_PI;
	vec3f dir(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
	dir = genFloat()>=0.5?dir:-dir;
	return dir;
	*/
	float z = 1.f - 2.f * genFloat();
	float r = sqrtf(max(0.f , 1.f - z * z));
	float phi = 2 * M_PI * genFloat();
	float x = r * cos(phi);
	float y = r * sin(phi);
	return vec3f(x , z , y);
}

vec3f RandGenerator::genHemisphericalDirection()
{
	float z = genFloat();
	float r = sqrtf(max(0.f , 1.f - z * z));
	float phi = 2 * M_PI * genFloat();
	float x = r * cos(phi);
	float y = r * sin(phi);
	return vec3f(x , z , y);
}

vec3f RandGenerator::genConcetricDisk()
{
	float phi , r;

	float a = 2 * genFloat() - 1;   /* (a,b) is now on [-1,1]^2 */
	float b = 2 * genFloat() - 1;

	if (a > -b)      /* region 1 or 2 */
	{
		if (a > b)   /* region 1, also |a| > |b| */
		{
			r = a;
			phi = (M_PI / 4.f) * (b / a);
		}
		else        /* region 2, also |b| > |a| */
		{
			r = b;
			phi = (M_PI / 4.f) * (2.f - (a / b));
		}
	}
	else            /* region 3 or 4 */
	{
		if(a < b)   /* region 3, also |a| >= |b|, a != 0 */
		{
			r = -a;
			phi = (M_PI / 4.f) * (4.f + (b / a));
		}
		else        /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
		{
			r = -b;

			if (b != 0)
				phi = (M_PI / 4.f) * (6.f - (a / b));
			else
				phi = 0;
		}
	}

	vec3f res;
	res.x = r * std::cos(phi);
	res.y = r * std::sin(phi);
	res.z = 0.f;
	return res;
}

vec3f RandGenerator::genHemiCosDirection(float expTerm, float *pdf)
{
	vec3f disk = vec3f(genFloat() , genFloat() , 0.f);

	float u1 = 2.f * M_PI * disk.x;
	float u2 = std::pow(disk.y , 1.f / (expTerm + 1.f));
	float u3 = std::sqrt(1.f - u2 * u2);

	vec3f res(std::cos(u1) * u3 , u2 , std::sin(u1) * u3);

	res.normalize();

	if (pdf)
		*pdf = (expTerm + 1.f) * std::pow(u2 , expTerm) * (0.5f / M_PI);

	return res;
	//printf("%.8f\n" , *pdf);
	
	/*
	float phi = acos(clampf(powf(genFloat(), 1/(expTerm+1)), -1, 1));
	float theta = genFloat()*2*M_PI;
	vec3f dir(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
	
	vec3f up = vec3f(0, 1, 0);
	vec3f axis = up.cross(normal);
	axis.normalize();
	if (axis.length() < 1e-6f)
		return res;

	float angle = acos(clampf(up.dot(normal), -1, 1));
	vec3f dir = vec3f(rotMat(axis, angle)*vec4<float>(res, 0));
	
	return dir;
	*/
}

/*
random_device RandGenerator::rd;
tr1::mt19937 RandGenerator::eng(19931004);
tr1::uniform_real_distribution<float> RandGenerator::dis_float(0, 1);
tr1::uniform_real_distribution<double> RandGenerator::dis_double(0, 1);

float RandGenerator::genFloat()
{
	return dis_float(eng);
}
*/