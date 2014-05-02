#pragma once

#include <vector>
#include <cmath>
#include "nvVector.h"
#include "UniformSphericalSampler.h"
using namespace nv;

typedef unsigned int uint;

struct CountQuery
{
	vec3f& pos;
	double count;

	CountQuery(vec3f& _pos) : pos(_pos) , count(0.f) {}
};

class CountHashGrid
{
public:
	~CountHashGrid()
	{
		clear();
	}

	void clear()
	{
		weights.clear();
		weights.shrink_to_fit();
	}

	void init(Scene& scene)
	{
		vec3f diag = scene.getDiagonal();
		mBBoxMin = scene.tree.root->boundingBox.minCoord;
		mBBoxMax = scene.tree.root->boundingBox.maxCoord;

		sizeNum = 100;

		mCellSize = max(diag[0] , max(diag[1] , diag[2])) / sizeNum;
		mInvCellSize = 1.f / mCellSize;
		cellArea = mCellSize * mCellSize;
		cellVolume = mCellSize * mCellSize * mCellSize;
		printf("cell size = %.8f\n" , mCellSize);
		
		totVolume = diag[0] * diag[1] * diag[2];

		weights.clear();
		weights.resize(sizeNum * sizeNum * sizeNum);
		memset(&weights[0] , 0 , sizeof(weights));

		sumContribs = 0;
	}

	template<typename tParticle>
	void addPhotons(const std::vector<tParticle> &aParticles , const int st , const int ed)
	{
		double energy;

		for (int i = 0; i < weights.size(); i++)
			weights[i] *= sumContribs;

		for(int i=st; i<ed; i++)
		{
			if (aParticles[i].ray->insideObject == NULL)
				continue;
			const vec3f &pos = aParticles[i].pos;
			vec3f totContrib = aParticles[i].dirContrib;
			int cellIndex = GetCellIndex(pos);
			if (cellIndex == -1)
			{
				printf("error index\n");
				continue;
			}
			
			energy = y(totContrib) * 1e-7f;
			sumContribs += energy;
			weights[cellIndex] += energy;
		}

		for (int i = 0; i < weights.size(); i++)
			weights[i] /= sumContribs;

		effectiveIndex.clear();
		effectiveWeights.clear();
		for (int i = 0; i < weights.size(); i++)
		{
			if (weights[i] < 1e-7f)
				continue;
			effectiveIndex.push_back(i);
			effectiveWeights.push_back(weights[i]);
			int N = effectiveWeights.size();
			effectiveWeights[N - 1] += effectiveWeights[N - 2];
		}
	}

	template<typename tQuery>
	void count(tQuery& aQuery)
	{
		const vec3f queryPos = aQuery.pos;
		const vec3f distMin = queryPos - mBBoxMin;
		const vec3f distMax = mBBoxMax - queryPos;
		for(int i=0; i<3; i++)
		{
			if(distMin[i] < 0.f) return;
			if(distMax[i] < 0.f) return;
		}

		const vec3f cellPt = mInvCellSize * distMin;
		const vec3f coordF(
			std::floor(cellPt.x),
			std::floor(cellPt.y),
			std::floor(cellPt.z));

		const int px = int(coordF.x);
		const int py = int(coordF.y);
		const int pz = int(coordF.z);

		const vec3f fractCoord = cellPt - coordF;

		const int pxo = px + (fractCoord.x < 0.5f ? -1 : +1);
		const int pyo = py + (fractCoord.y < 0.5f ? -1 : +1);
		const int pzo = pz + (fractCoord.z < 0.5f ? -1 : +1);

		for(int j=0; j<7; j++)
		{
			int cellIndex;
			switch(j)
			{
				case 0: cellIndex = (GetCellIndex(vec3i(px , py , pz ))); break;
				case 1: cellIndex = (GetCellIndex(vec3i(px , py , pzo))); break;
				case 2: cellIndex = (GetCellIndex(vec3i(px , pyo, pz ))); break;
				case 3: cellIndex = (GetCellIndex(vec3i(px , pyo, pzo))); break;
				case 4: cellIndex = (GetCellIndex(vec3i(pxo, py , pz ))); break;
				case 5: cellIndex = (GetCellIndex(vec3i(pxo, py , pzo))); break;
				case 6: cellIndex = (GetCellIndex(vec3i(pxo, pyo, pz ))); break;
				case 7: cellIndex = (GetCellIndex(vec3i(pxo, pyo, pzo))); break;
			}
			if (cellIndex == -1)
				continue;
			aQuery.count += weights[cellIndex] / 8;
		}
	}

	void print(FILE* fp)
	{
		fprintf(fp , "============ one iter, sum = %.8f ============\n" , sumContribs);
		for (int i = 0; i < weights.size(); i++)
		{
			if (weights[i] <= 0)
				continue;
			fprintf(fp , "i = %d, weight = %.8f\n" , i , weights[i]);
		}
	}

public:
	int GetCellIndex(const vec3i &aCoord) const
	{
		uint x = uint(aCoord.x);
		uint y = uint(aCoord.y);
		uint z = uint(aCoord.z);

		if (x < 0 || x >= sizeNum || y < 0 || y >= sizeNum || z < 0 || z >= sizeNum)
		{
			//printf("error: x = %d, y = %d, z = %d\n" , x , y , z);
			return -1;
		}

		return int(((x * sizeNum * sizeNum) + (y * sizeNum) + z));
	}

	int GetCellIndex(const vec3f &aPoint) const
	{
		const vec3f distMin = aPoint - mBBoxMin;

		vec3f coordF(
			std::floor(mInvCellSize * distMin.x),
			std::floor(mInvCellSize * distMin.y),
			std::floor(mInvCellSize * distMin.z));

		for (int i = 0; i < 3; i++)
			coordF[i] = clamp(coordF[i] , 0 , sizeNum - 1);

		const vec3i coordI  = vec3i(int(coordF.x), int(coordF.y), int(coordF.z));

		return GetCellIndex(coordI);
	}
	
	vec3f cellIndexToPosition(const int &index)
	{
		int x = index / (sizeNum * sizeNum);
		int y = index / sizeNum % sizeNum;
		int z = index % sizeNum;

		vec3f corner = mBBoxMin + vec3f(mCellSize , 0 , 0) * x + 
			vec3f(0 , mCellSize , 0) * y + vec3f(0 , 0 , mCellSize) * z;

		vec3f res = corner + vec3f(mCellSize , 0 , 0) * RandGenerator::genFloat() +
			vec3f(0 , mCellSize , 0) * RandGenerator::genFloat() + vec3f(0 , 0 , mCellSize) * RandGenerator::genFloat();

		return res;
	}

	vec3f getRandomPosition(float &pdf)
	{
		float rnd = RandGenerator::genFloat();
		unsigned index = (lower_bound(effectiveWeights.begin(), effectiveWeights.end(), rnd)-effectiveWeights.begin());
		if(index >= effectiveWeights.size() - 1)
		{
			index = effectiveWeights.size() - 1;
			pdf = 1.f - effectiveWeights[index];
		}
		else
		{
			pdf = effectiveWeights[index + 1] - effectiveWeights[index];
		}
		return cellIndexToPosition(effectiveIndex[index]);
	}

	Ray volumeEmit(Scene *scene)
	{
		Ray ray;
		ray.contactObject = NULL;
		ray.contactObjectTriangleID = -1;
		ray.origin = getRandomPosition(ray.originProb);

		UniformSphericalSampler uniformSphericalSampler;

		LocalFrame lf;

		ray.direction = uniformSphericalSampler.genSample(lf);

		ray.insideObject = scene->findInsideObject(ray, ray.contactObject);

		ray.current_tid = scene->getContactTreeTid(ray);
		ray.color = vec3f(1, 1, 1);

		ray.directionProb = uniformSphericalSampler.getProbDensity(lf , ray.direction);
		
		ray.directionSampleType = ray.originSampleType = Ray::RANDOM;

		if(!scene->usingGPU())
		{
			Scene::ObjSourceInformation osi;
			NoSelfIntersectionCondition condition(scene, ray);
			float dist = scene->intersect(ray, osi, &condition);
			if(dist > 0)
			{
				ray.intersectDist = dist;
				ray.intersectObject = scene->objects[osi.objID];
				ray.intersectObjectTriangleID = osi.triangleID;
			}
		}
		return ray;
	}

public:
	vec3f mBBoxMin;
	vec3f mBBoxMax;
	std::vector<double> weights;

	std::vector<double> effectiveWeights;
	std::vector<int> effectiveIndex;

	float mCellSize;
	float mInvCellSize;

	double sumContribs;
	double cellArea;
	double cellVolume;
	double totVolume;
	int sizeNum;
};

