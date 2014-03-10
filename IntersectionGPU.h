#pragma once

#include "SimpleShape.h"
#include "macros.h"
#include "KDTree.h"
#include "Shader.h"
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#include <GL/glew.h>
#include <sstream>

class IntersectionGPU
{
public:
	struct Ray
	{
		vec3f origin;
		vec3f direction;
		int last_tid;
		Ray(){}
		Ray(const vec3f& origin, const vec3f& direction)
		{
			this->origin = origin;
			this->direction = direction;
		}
	};
protected:
	static int maxTexSize;
	static bool renderToFBO;
	static int width, height;
	static GLuint fboID, color_rboID, depth_rboID;
	
	static Shader shader;
	static void initFBO(int w = 0, int h = 0);
	static void cleanUpFBO();
	static void render();
	static vector<vec4f> queryOneBatch(const vector<Ray>& rayList);
	
public:
	static vector<SimpleShape*> shapes;
	static vector<Ray> rayList;
	static void loadShapesToShader();
	static void loadKDTreeToShader(const KDTree& tree);
	static void init();
	static vector<vec4f> query(const vector<Ray>& rayList);
};