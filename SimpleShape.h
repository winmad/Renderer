#pragma once
#include <string>
#include <vector>
#include "RandGenerator.h"
#include "smallFuncs.h"
#include "LocalFrame.h"

using namespace std;
using namespace nv;



#define BUFFERSIZE 1024

#define min(a,b) a<b?a:b
#define max(a,b) a>b?a:b

class SimpleShape
{
public:
	string name;

	vec3f minCoord, maxCoord;
	//v ID
	vector<vec3f> vertexList;

	//vn ID
	vector<vec3f> vertexNormalList;

	//vt ID
	vector<vec3f> vertexTexCoordList;

	//face ID
	vector<vec3ui> faceVertexIndexList;
	vector<vec3ui> faceVertexNormalIndexList;
	vector<vec3ui> faceVertexTexCoordIndexList;

	matrix4<float> transform;

	void setTransform(const matrix4<float>& mat) { transform = mat; }
	matrix4<float> getTransform() const { return transform; }
	SimpleShape(){};
	SimpleShape(const string &fileName, bool normalize = true){ loadShape(fileName, normalize); }
	vec3f getVertexPosition(int vi) const { return vertexList[vi]; }
	vec3ui getVertexIndices(int ti) const { return faceVertexIndexList[ti]; }
	vec3f getWorldVertexPosition(int vi) const { return vec3f(transform*vec4<float>(vertexList[vi], 1)); }
	virtual vec3f getTexCoord(unsigned fi, const vec3f& position) const;
	virtual vec3f getWorldNormal(unsigned fi, const vec3f& position, bool flat = false) const;
	virtual LocalFrame getAutoGenWorldLocalFrame(unsigned fi, const vec3f& position, bool flat = false) const;
	unsigned getVertexNum() const { return vertexList.size(); }
	unsigned getTriangleNum() const { return faceVertexIndexList.size(); }
	float getTriangleArea(unsigned ti) const;
	vec3f genRandTrianglePosition(unsigned ti) const;
	vec3f getCenter() const;
	float getDiagLen() const;
	void getBoundingBox(vec3f &minCoord, vec3f &maxCoord);
	void loadShape(const string &fileName, bool normalize = true, vector<SimpleShape*>* ss = NULL);
	matrix4<float> unitize();
	void saveShape(const string &fileName);
	void applyTransform();
};