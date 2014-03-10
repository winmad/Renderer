#pragma once
#include "nvVector.h"
#include "nvMatrix.h"
#include "EXT/rapidxml/rapidxml.hpp"

using namespace rapidxml;
using namespace nv;

matrix4<float> readMatrix(const char* value);

vec3f readVec(const char* value);