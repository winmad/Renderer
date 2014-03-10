#include "nvVector.h"

using namespace nv;

#define EPSILON (1e-5)
#define EPS_RAY (1e-5)
#define COS_TERM_MIN (1e-4)

#define M_PI 3.14159265358979
#define max2(a, b) a > b ? a : b
#define min2(a, b) a < b ? a : b

typedef vec2<double> vec2d;
typedef vec3<double> vec3d;
