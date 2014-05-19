#pragma once

#include <vector>
#include <cmath>
#include "nvVector.h"
#include "UniformSphericalSampler.h"
#include "HGPhaseSampler.h"
#include "LocalFrame.h"
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

	void preprocess(Scene& scene)
	{
		effectiveIndex.clear();
		effectiveWeights.clear();

		Ray ray;
		float sumWeights = 0.f;

		omp_lock_t lock;
		omp_init_lock(&lock);

#pragma omp parallel for
		for (int i = 0; i < sizeNum * sizeNum * sizeNum; i++)
		{
			bool isInside = false;
			for (int j = 0; j < 8; j++)
			{
				vec3f offset;
				offset.x = (j & 1);
				offset.y = ((j >> 1) & 1);
				offset.z = ((j >> 2) & 1);
				ray.origin = cellIndexToPosition(i , offset);
				ray.direction = RandGenerator::genSphericalDirection();
				SceneObject *insideObject = scene.findInsideObject(ray);
				if (insideObject && insideObject->isVolumeric())
				{
					isInside = true;
					break;
				}
			}
			if (!isInside)
			{
				vec3f offset;
				offset.x = offset.y = offset.z = 0.5f;
				ray.origin = cellIndexToPosition(i , offset);
				ray.direction = RandGenerator::genSphericalDirection();
				SceneObject *insideObject = scene.findInsideObject(ray);
				if (insideObject && insideObject->isVolumeric())
					isInside = true;
			}

			if (isInside)
			{
				omp_set_lock(&lock);
				sumWeights += 1.f;
				effectiveIndex.push_back(i);
				effectiveWeights.push_back(sumWeights);
				omp_unset_lock(&lock);
			}
		}

		for (int i = 0; i < effectiveWeights.size(); i++)
			effectiveWeights[i] /= sumWeights;

		totVolume = sumWeights * mCellSize * mCellSize * mCellSize;

		omp_destroy_lock(&lock);
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
			vec3f totContrib = aParticles[i].throughput;
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
			//effectiveWeights.push_back(1.f);
		}

		//for (int i = 0; i < effectiveWeights.size(); i++)
		//	effectiveWeights[i] /= (Real)effectiveWeights.size();
		for (int i = 1; i < effectiveWeights.size(); i++)
			effectiveWeights[i] += effectiveWeights[i - 1];
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

		int N = 1;
		for(int j=0; j<N; j++)
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
			aQuery.count += weights[cellIndex] / (float)N;
		}
	}

	void print(FILE* fp)
	{
		fprintf(fp , "============ one iter============\n");
		for (int i = 0; i < effectiveWeights.size(); i++)
		{
			fprintf(fp , "index = %d, accuWeight = %.8f\n" , effectiveIndex[i] , effectiveWeights[i]);
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
	
	// offset in [(0,0,0) , (1,1,1)]
	vec3f cellIndexToPosition(const int &index , const vec3f& offset)
	{
		int x = index / (sizeNum * sizeNum);
		int y = index / sizeNum % sizeNum;
		int z = index % sizeNum;

		vec3f corner = mBBoxMin + vec3f(mCellSize , 0 , 0) * x + 
			vec3f(0 , mCellSize , 0) * y + vec3f(0 , 0 , mCellSize) * z;

		vec3f res = corner + vec3f(mCellSize , mCellSize , mCellSize) * offset;

		return res;
	}

	vec3f getRandomPosition(float &pdf)
	{
		float rnd = RandGenerator::genFloat();
		unsigned index = (lower_bound(effectiveWeights.begin(), effectiveWeights.end(), rnd)-effectiveWeights.begin());
		if(index == 0)
		{
			pdf = effectiveWeights[index];
		}
		else
		{
			pdf = effectiveWeights[index] - effectiveWeights[index - 1];
		}
		vec3f rndOffset = vec3f(RandGenerator::genFloat() , RandGenerator::genFloat() , RandGenerator::genFloat());
		vec3f res = cellIndexToPosition(effectiveIndex[index] , rndOffset);
		return res;
	}

	Ray volumeEmit(Scene *scene)
	{
		for (;;)
		{
		Ray ray;
		ray.contactObject = NULL;
		ray.contactObjectTriangleID = -1;
		ray.origin = getRandomPosition(ray.originProb);

		// not sure
		ray.originProb *= mInvCellSize * mInvCellSize * mInvCellSize;
		//ray.originProb = 1.f / scene->getTotalVolume();

		//printf("%.8f , %.8f\n" , ray.originProb , 1.f / scene->getTotalVolume());

		ray.direction = RandGenerator::genSphericalDirection();
		ray.insideObject = scene->findInsideObject(ray, ray.contactObject);
		if (ray.insideObject == NULL)
		{
			continue;
		}

		HGPhaseSampler hgSampler(ray.insideObject->getG());
		LocalFrame lf;
		lf.buildFromNormal(RandGenerator::genSphericalDirection());
		/*
		ray.direction = RandGenerator::genSphericalDirection();
		ray.directionProb = 0.25f / M_PI;
		*/
		ray.direction = hgSampler.genSample(lf);
		ray.directionProb = hgSampler.getProbDensity(lf , ray.direction);
		//printf("%.8f\n" , ray.directionProb);

		ray.current_tid = scene->getContactTreeTid(ray);
		ray.color = vec3f(1, 1, 1);

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

