#pragma once

#include "nvVector.h"
#include <vector>
#include <limits>

typedef float Real;
typedef nv::vec3f Vector3;

#define INF std::numeric_limits<float>::infinity()

struct KdNode
{
	Real splitPos;
	uint32_t splitAxis : 2;
	uint32_t hasLeftChild : 1 , rightChild : 29;

	void init(Real pos , uint32_t a)
	{
		splitPos = pos;
		splitAxis = a;
		rightChild = (1 << 29) - 1;
		hasLeftChild = 0;
	}

	void initLeaf()
	{
		splitAxis = -1;
		rightChild = (1 << 29) - 1;
		hasLeftChild = 0;
	}
};

template<typename NodeData>
class PointKDTree
{
public:
	KdNode *nodes;
	NodeData *nodeData;
	uint32_t totNodes , nextFreeNode;

	PointKDTree(const std::vector<NodeData>& data);

	~PointKDTree()
	{
		delete[] nodes;
		delete[] nodeData;
	}

	void buildTree(uint32_t nodeNum , int st , int ed ,
		const NodeData **buildNodes);

	template<typename Query>
	void searchInRadius(uint32_t nodeNum , const Vector3& pos ,
		Real radius , Query& query);

	template<typename Query>
	void searchKnn(uint32_t nodeNum , const Vector3& pos , 
		Query& query);
};

template<typename NodeData>
PointKDTree<NodeData>::PointKDTree(const std::vector<NodeData>& data)
{
	totNodes = data.size();
	nextFreeNode = 1;
	nodes = new KdNode[totNodes];
	nodeData = new NodeData[totNodes];
	std::vector<const NodeData*> buildNodes(totNodes , NULL);
	for (uint32_t i = 0; i < totNodes; i++)
	{
		buildNodes[i] = &data[i];
	}
	buildTree(0 , 0 , totNodes , &buildNodes[0]);
}

template<typename NodeData>
struct CompareNode
{
	int axis;

	CompareNode(int axis) : axis(axis) {}

	bool operator()(const NodeData *d1 , const NodeData *d2) const
	{
		return (d1->pos[axis] == d2->pos[axis]) ? (d1 < d2) : 
			(d1->pos[axis] < d2->pos[axis]);
	}
};

template<typename NodeData>
void PointKDTree<NodeData>::buildTree(uint32_t nodeNum , 
	int st , int ed , const NodeData **buildNodes)
{
	if (st + 1 == ed)
	{
		nodes[nodeNum].initLeaf();
		nodeData[nodeNum] = *buildNodes[st];
		return;
	}

	Vector3 l(INF) , r(-INF);
	for (int i = st; i < ed; i++)
	{
		l.x = std::min(l.x , buildNodes[i]->pos.x);
		l.y = std::min(l.y , buildNodes[i]->pos.y);
		l.z = std::min(l.z , buildNodes[i]->pos.z);
		r.x = std::max(r.x , buildNodes[i]->pos.x);
		r.y = std::max(r.y , buildNodes[i]->pos.y);
		r.z = std::max(r.z , buildNodes[i]->pos.z);
	}

	Vector3 diag = r - l;
	Real tmp = -INF;
	int splitAxis = -1;
	for (int i = 0; i <= 2; i++)
	{
		if (tmp < diag[i])
		{
			tmp = diag[i];
			splitAxis = i;
		}
	}
	int splitPos = (st + ed) / 2;
	std::nth_element(&buildNodes[st] , &buildNodes[splitPos] ,
		&buildNodes[ed] , CompareNode<NodeData>(splitAxis));

	nodes[nodeNum].init(buildNodes[splitPos]->pos[splitAxis] , splitAxis);
	nodeData[nodeNum] = *buildNodes[splitPos];

	if (st < splitPos)
	{
		nodes[nodeNum].hasLeftChild = 1;
		uint32_t childNum = nextFreeNode++;
		buildTree(childNum , st , splitPos , buildNodes);
	}
	if (splitPos + 1 < ed)
	{
		nodes[nodeNum].rightChild = nextFreeNode++;
		buildTree(nodes[nodeNum].rightChild , splitPos + 1 , ed , buildNodes);
	}
}

template<typename NodeData> 
template<typename Query>
void PointKDTree<NodeData>::searchInRadius(uint32_t nodeNum , const Vector3& pos ,
	Real radius , Query& query)
{
	KdNode *node = &nodes[nodeNum];

	int axis = node->splitAxis;

	Real delta = std::abs(pos[axis] - node->splitPos);

	if (axis != -1)
	{
		if (pos[axis] <= node->splitPos)
		{
			if (node->hasLeftChild)
				searchInRadius(nodeNum + 1 , pos , radius , query);
			if (delta < radius && node->rightChild < totNodes)
				searchInRadius(node->rightChild , pos , radius , query);
		}
		else
		{
			if (node->rightChild < totNodes)
				searchInRadius(node->rightChild , pos , radius , query);
			if (delta < radius && node->hasLeftChild)
				searchInRadius(nodeNum + 1 , pos , radius , query);
		}
	}

	Vector3 d = pos - nodeData[nodeNum].pos;
	Real dis = d.length();
	if (dis < radius)
		query.process(nodeData[nodeNum]);
}

template<typename NodeData>
template<typename Query>
void PointKDTree<NodeData>::searchKnn(uint32_t nodeNum , const Vector3& pos ,
	Query& query)
{
	KdNode *node = &nodes[nodeNum];

	int axis = node->splitAxis;

	Real delta2 = powf(pos[axis] - node->splitPos, 2);

	if (axis != -1)
	{
		if (pos[axis] <= node->splitPos)
		{
			if (node->hasLeftChild)
				searchKnn(nodeNum + 1 , pos , query);
			if (delta2 < query.maxSqrDis && node->rightChild < totNodes)
				searchKnn(node->rightChild , pos , query);
		}
		else
		{
			if (node->rightChild < totNodes)
				searchKnn(node->rightChild , pos , query);
			if (delta2 < query.maxSqrDis && node->hasLeftChild)
				searchKnn(nodeNum + 1 , pos , query);
		}
	}

	Vector3 d = pos - nodeData[nodeNum].pos;
	Real dist2 = powf(d.length(), 2);
	if (dist2 < query.maxSqrDis)
		query.process(&nodeData[nodeNum] , dist2);
}
