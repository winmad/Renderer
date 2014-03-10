#pragma once
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "nvVector.h"

using namespace std;
using namespace nv;



class KDTree
{
public:
	struct Ray
	{
		vec3f origin;
		vec3f direction;
	};

	struct Triangle
	{
		vec3ui vertexIndices;
		void* sourceInformation;
	};
	enum Strategy{ LOOP, BEST } strategy;

	class Condition
	{
	public:
		virtual bool legal(const Ray& ray, const Triangle& tri, const float dist) const
		{
			return true;
		}
	};

private:
	struct Node;
	Node* root;
	unsigned maxLeafTriangleNum;
	unsigned maxDepth;
	void splitLoop(Node* node, unsigned dim, unsigned depth = 0);
	void splitBest(Node* node, unsigned depth = 0);
	void destroy(Node* node);
	float intersect(const Ray& ray, const Node* node, unsigned& triangleID, const Condition* condition) const;
	float intersect(const Ray& ray, const Triangle& tri, const Condition* condition) const;
	void serializeForGPU(Node* node, int parent, vector<vec4f>& nodes, vector<vec4f>& nodes_minCoords, vector<vec4f>& nodes_maxCoords, vector<vec4f>& leaf_v1, vector<vec4f>& leaf_v2, vector<vec4f>& leaf_v3) const;
public:
	vector<vec3f> vertexPositionList;
	vector<Triangle> triangleList;

	float intersect(const Ray& ray, unsigned& triangleID, const Condition* condition = NULL) const;

	KDTree();
	void build(Strategy strategy = BEST);
	void destroy();
	void serializeForGPU(vector<vec4f>& nodes, vector<vec4f>& nodes_minCoords, vector<vec4f>& nodes_maxCoords, vector<vec4f>& leaf_v1, vector<vec4f>& leaf_v2, vector<vec4f>& leaf_v3) const;
	~KDTree();

};

