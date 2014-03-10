#pragma once
#include "GL/glew.h"
#include "textfile.h"
#include "nvVector.h"
#include "nvMatrix.h"
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;
using namespace nv;

class Shader
{
protected:
	vector<GLuint> texIDs;
	unsigned maxTexUnit;
	unordered_map<string, GLuint> progName2progID;
	string currentProgName;
public:
	void init();
	bool createProgram(const string& progName, const vector<string>& vertFileNameList, const vector<string>& fragFileNameList);
	bool createProgram(const string& progName, const string& vertFileName, const string& fragFileName)
	{
		return createProgram(progName, vector<string>(1, vertFileName), vector<string>(1, fragFileName));
	}
	bool useProgram(const string& progName);
	bool deleteProgram(const string& progName);
	void setUniform(const string& progName, const string& varName, const int v);
	void setUniform(const string& progName, const string& varName, const float v);
	void setUniform(const string& progName, const string& varName, const vec2f& v);
	void setUniform(const string& progName, const string& varName, const vec3f& v);
	void setUniform(const string& progName, const string& varName, const vec4f& v);
	void setUniform(const string& progName, const string& varName, const matrix4<float>& v);
	void setUniform(const string& varName, const int v) { setUniform(currentProgName, varName, v); }
	void setUniform(const string& varName, const float v) { setUniform(currentProgName, varName, v); }
	void setUniform(const string& varName, const vec2f& v) { setUniform(currentProgName, varName, v); }
	void setUniform(const string& varName, const vec3f& v) { setUniform(currentProgName, varName, v); }
	void setUniform(const string& varName, const vec4f& v) { setUniform(currentProgName, varName, v); }
	void setUniform(const string& varName, const matrix4<float>& v) { setUniform(currentProgName, varName, v); }
	void setTexArray1D(const string& varName, const vector<vec4f>& tex);
	void setTexArray2D(const string& varName, const vector<vec4f>& tex, int width = 1);
	void setTexArray2D(const string& varName, const vector<vector<vec4f>>& tex);
	void deleteLastTexture()
	{
		if(texIDs.size())
		{
			glDeleteTextures(1, &texIDs.back());
			texIDs.pop_back();
			maxTexUnit--;
		}
	}
	void clearTextures()
	{
		glDeleteTextures(texIDs.size(), texIDs.data());
		texIDs.clear();
		maxTexUnit = 0;
	}
	~Shader(void);
};

