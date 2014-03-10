#include "StdAfx.h"
#include "KDTree.h"

template<bool min> vec3f boundvec(const vec3f& v1, const vec3f& v2)
{
	vec3f result;
	if(min)
	{
		result.x = v1.x < v2.x ? v1.x : v2.x;
		result.y = v1.y < v2.y ? v1.y : v2.y;
		result.z = v1.z < v2.z ? v1.z : v2.z;
	}
	else
	{
		result.x = v1.x > v2.x ? v1.x : v2.x;
		result.y = v1.y > v2.y ? v1.y : v2.y;
		result.z = v1.z > v2.z ? v1.z : v2.z;
	}
	return result;
}

KDTree::KDTree()
{
	maxLeafTriangleNum = 1;
	maxDepth = -1;
	root = NULL;
}

struct KDTree::Node
{
	Node* left;
	Node* right;
	vector<unsigned> triangleIndices;
	struct BoundingBox
	{
		vec3f minCoord;
		vec3f maxCoord;
		float intersect(const KDTree::Ray& ray) const;
	} boundingBox;
};



float KDTree::Node::BoundingBox::intersect(const Ray& ray) const
{
	vec3f intPoint;
	float minDist = -1;
	for(unsigned i=0; i<6; i++)
	{
		if(ray.direction[i%3] == 0)
			continue;
		float dist = (((float*)this)[i] - ray.origin[i%3]) / ray.direction[i%3];
		intPoint = ray.origin + dist*ray.direction;
		bool inter = dist >= 0;
		inter &= minCoord[(i+1)%3] <= intPoint[(i+1)%3];
		inter &= maxCoord[(i+1)%3] >= intPoint[(i+1)%3];
		inter &= minCoord[(i+2)%3] <= intPoint[(i+2)%3];
		inter &= maxCoord[(i+2)%3] >= intPoint[(i+2)%3];
		minDist = (inter && (dist < minDist || minDist < 0)) ? dist : minDist;
	}
	return minDist;
}

void KDTree::splitLoop(Node* node, unsigned dim, unsigned depth)
{
	node->left = node->right = NULL;
	if(node->triangleIndices.size() == 0)
		return;
	node->boundingBox.minCoord = node->boundingBox.maxCoord = 
		vertexPositionList[triangleList[node->triangleIndices[0]].vertexIndices[0]];

	unordered_set<unsigned> nodeVertexIndices;

	for(unsigned i=0; i<node->triangleIndices.size(); i++)
	{
		for(unsigned k=0; k<3; k++)
		{
			vec3f &p = vertexPositionList[triangleList[node->triangleIndices[i]].vertexIndices[k]];
			node->boundingBox.minCoord = boundvec<1>(p, node->boundingBox.minCoord);
			node->boundingBox.maxCoord = boundvec<0>(p, node->boundingBox.maxCoord);
			nodeVertexIndices.insert(triangleList[node->triangleIndices[i]].vertexIndices[k]);
		}
	}

	if(depth > maxDepth)
		return;

	if(node->triangleIndices.size() <= maxLeafTriangleNum)
		return;

	vector<float> coordList;
	for(unordered_set<unsigned>::iterator it=nodeVertexIndices.begin(); it!=nodeVertexIndices.end(); it++)
	{
		coordList.push_back(vertexPositionList[*it][dim]);
	}

	nth_element(coordList.begin(), coordList.begin() + (coordList.size()>>1), coordList.end());

	nodeVertexIndices.clear();

	float splitValue = coordList[coordList.size()>>1];

	node->left = new Node;
	node->right = new Node;

	unsigned leftCount = 0;

	for(unsigned i=0; i<node->triangleIndices.size(); i++)
	{
		unsigned mask = 0;
		for(unsigned k=0; k<3; k++)
		{
			vec3f &p = vertexPositionList[triangleList[node->triangleIndices[i]].vertexIndices[k]];
			mask += (p[dim] <= splitValue) << k;
		}
		if(mask == 0)
		{
			node->right->triangleIndices.push_back(node->triangleIndices[i]);
		}
		else
		{
			node->left->triangleIndices.push_back(node->triangleIndices[i]);
			leftCount++;
		}
	}

	if(leftCount == 0 || leftCount == node->triangleIndices.size())
	{
		delete node->left;
		delete node->right;
		node->left = node->right = NULL;
		return;
	}

	splitLoop(node->left, (dim+1)%3, depth + 1);
	splitLoop(node->right, (dim+1)%3, depth + 1);
}

void KDTree::splitBest(Node* node, unsigned depth)
{
	float bestRatio = 10;
	unsigned bestDim;
	float bestSplitValue;
	for(unsigned dim=0; dim<3; dim++)
	{
		node->left = node->right = NULL;
		if(node->triangleIndices.size() == 0)
			return;
		node->boundingBox.minCoord = node->boundingBox.maxCoord = 
			vertexPositionList[triangleList[node->triangleIndices[0]].vertexIndices[0]];

		unordered_set<unsigned> nodeVertexIndices;

		for(unsigned i=0; i<node->triangleIndices.size(); i++)
		{
			for(unsigned k=0; k<3; k++)
			{
				vec3f &p = vertexPositionList[triangleList[node->triangleIndices[i]].vertexIndices[k]];
				node->boundingBox.minCoord = boundvec<1>(p, node->boundingBox.minCoord);
				node->boundingBox.maxCoord = boundvec<0>(p, node->boundingBox.maxCoord);
				nodeVertexIndices.insert(triangleList[node->triangleIndices[i]].vertexIndices[k]);
			}
		}

		if(depth > maxDepth)
			return;

		if(node->triangleIndices.size() <= maxLeafTriangleNum)
			return;

		vector<float> coordList;
		for(unordered_set<unsigned>::iterator it=nodeVertexIndices.begin(); it!=nodeVertexIndices.end(); it++)
		{
			coordList.push_back(vertexPositionList[*it][dim]);
		}

		nth_element(coordList.begin(), coordList.begin() + (coordList.size()>>1), coordList.end());

		nodeVertexIndices.clear();

		float splitValue = coordList[coordList.size()>>1];

		unsigned leftCount = 0;

		for(unsigned i=0; i<node->triangleIndices.size(); i++)
		{
			unsigned mask = 0;
			
			for(unsigned k=0; k<3; k++)
			{
				vec3f &p = vertexPositionList[triangleList[node->triangleIndices[i]].vertexIndices[k]];
				mask += (p[dim] <= splitValue) << k;
			}

			if(mask != 0)
			{
				leftCount++;
			}
		}
		float ratio = leftCount/float(node->triangleIndices.size());
		if(abs(ratio-0.5) < abs(bestRatio - 0.5))
		{
			bestRatio = ratio;
			bestDim = dim;
			bestSplitValue = splitValue;
		}
	}

	if(bestRatio == 0 || bestRatio == 1)
	{
		node->left = node->right = NULL;
		return;
	}

	node->left = new Node;
	node->right = new Node;

	for(unsigned i=0; i<node->triangleIndices.size(); i++)
	{
		unsigned mask = 0;
		for(unsigned k=0; k<3; k++)
		{
			vec3f &p = vertexPositionList[triangleList[node->triangleIndices[i]].vertexIndices[k]];
			mask += (p[bestDim] <= bestSplitValue) << k;
		}
		if(mask == 0)
		{
			node->right->triangleIndices.push_back(node->triangleIndices[i]);
		}
		else
		{
			node->left->triangleIndices.push_back(node->triangleIndices[i]);
		}
	}

	splitBest(node->left, depth + 1);
	splitBest(node->right, depth + 1);
}

void KDTree::build(Strategy strategy)
{
	destroy(root);
	root = new Node;
	for(unsigned i=0; i<triangleList.size(); i++)
		root->triangleIndices.push_back(i);

	switch(strategy)
	{
	case LOOP:
		splitLoop(root, 0);
		break;
	case BEST:
		splitBest(root);
		break;
	}
}

void KDTree::destroy(Node* node)
{
	if(node == NULL)
		return;
	if(node->left)
		destroy(node->left);
	if(node->right)
		destroy(node->right);
	delete node;
}

void KDTree::destroy()
{
	triangleList.clear();
	vertexPositionList.clear();
	destroy(root);
	root = NULL;
}

float KDTree::intersect(const Ray& ray, const Triangle& tri, const Condition* condition) const
{
	vec3f vs[3], lines[3], r2v, n;
	for(unsigned i=0; i<3; i++)
		vs[i] = vertexPositionList[tri.vertexIndices[i]];

	unsigned mask_pos = 0;
	unsigned mask_neg = 0;
	for(unsigned i=0; i<3; i++)
	{
		lines[i] = vs[(i+1)%3] - vs[i];
		r2v = vs[i] - ray.origin;
		n = lines[i].cross(r2v);
		mask_pos += (ray.direction.dot(n) >= 0) << i;
		mask_neg += (ray.direction.dot(n) <= 0) << i;
	}
	n = lines[0].cross(lines[1]);
	n.normalize();
	float orthoDist = n.dot(r2v);
	float d = n.dot(ray.direction);
	if(d == 0 || orthoDist/d < 0 || (mask_pos!=7 && mask_neg != 7))
		return -1;
	float dist = orthoDist/d;
	if(condition && !condition->legal(ray, tri, dist))
		return -1;
	return dist;
}

float KDTree::intersect(const Ray& ray, const Node* node, unsigned& triangleID, const Condition* condition) const
{
	if(node->boundingBox.intersect(ray) < 0)
		return -1;
	if(node->left == NULL && node->right == NULL)
	{
		float minDist = -1;
		for(unsigned ti=0; ti<node->triangleIndices.size(); ti++)
		{
			float dist = intersect(ray, triangleList[node->triangleIndices[ti]], condition);
			if(dist >= 0 && (dist < minDist || minDist < 0))
			{
				minDist = dist;
				triangleID = node->triangleIndices[ti];
			}
		}
		return minDist;
	}
	float ld, rd;
	unsigned ltid = 0, rtid = 0;
	ld = intersect(ray, node->left, ltid, condition);
	rd = intersect(ray, node->right, rtid, condition);
	if(ld >= 0 && rd >= 0)
	{
		triangleID = ld < rd ? ltid : rtid;
		return ld < rd ? ld : rd;
	}
	else if(rd < 0)
	{
		triangleID = ltid;
		return ld;
	}
	else if(ld < 0)
	{
		triangleID = rtid;
		return rd;
	}
	return -1;
}

float KDTree::intersect(const Ray& ray, unsigned& triangleID, const Condition* condition) const
{
	float dist = intersect(ray, root, triangleID, condition);
	return dist;
}

KDTree::~KDTree()
{
	destroy();
}

void KDTree::serializeForGPU(Node* node, int parent, vector<vec4f>& nodes, vector<vec4f>& nodes_minCoords, vector<vec4f>& nodes_maxCoords, vector<vec4f>& leaf_v1, vector<vec4f>& leaf_v2, vector<vec4f>& leaf_v3) const
{
	nodes_minCoords.push_back(vec4f(node->boundingBox.minCoord.x, node->boundingBox.minCoord.y, node->boundingBox.minCoord.z, 0));
	nodes_maxCoords.push_back(vec4f(node->boundingBox.maxCoord.x, node->boundingBox.maxCoord.y, node->boundingBox.maxCoord.z, 0));


	// nodes (parentIndex, rightChildIndex, leafStartIndex, leafEndIndex)
	if(node->left && node->right)
	{
		nodes.push_back(vec4f(parent, -1, -1, -1));
		int newParent = nodes.size()-1;
		serializeForGPU(node->left, newParent, nodes, nodes_minCoords, nodes_maxCoords, leaf_v1, leaf_v2, leaf_v3);
		nodes[newParent].y = nodes.size(); // fill right child
		serializeForGPU(node->right, newParent, nodes, nodes_minCoords, nodes_maxCoords, leaf_v1, leaf_v2, leaf_v3);
	}
	else // in leaf
	{
		int si = leaf_v1.size();
		for(unsigned i=0; i<node->triangleIndices.size(); i++)
		{
			const Triangle& tri = triangleList[node->triangleIndices[i]];
			vec3f v1 = vertexPositionList[tri.vertexIndices.x];
			vec3f v2 = vertexPositionList[tri.vertexIndices.y];
			vec3f v3 = vertexPositionList[tri.vertexIndices.z];
			leaf_v1.push_back(vec4f(v1.x, v1.y, v1.z, node->triangleIndices[i]));
			leaf_v2.push_back(vec4f(v2.x, v2.y, v2.z, 0));
			leaf_v3.push_back(vec4f(v3.x, v3.y, v3.z, 0));
		}
		nodes.push_back(vec4f(parent, -1, si, leaf_v1.size()));
	}
	
}

void KDTree::serializeForGPU(vector<vec4f>& nodes, vector<vec4f>& nodes_minCoords, vector<vec4f>& nodes_maxCoords, vector<vec4f>& leaf_v1, vector<vec4f>& leaf_v2, vector<vec4f>& leaf_v3) const
{
	// nodes (parentIndex, rightChildIndex, leafStartIndex, leafEndIndex)
	serializeForGPU(root, -1, nodes, nodes_minCoords, nodes_maxCoords, leaf_v1, leaf_v2, leaf_v3);
}