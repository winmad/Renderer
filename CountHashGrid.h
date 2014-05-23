#pragma once

#include <vector>
#include <cmath>
#include <omp.h>
#include "nvVector.h"
#include "UniformSphericalSampler.h"
#include "HGPhaseSampler.h"
#include "LocalFrame.h"
using namespace nv;

typedef unsigned int uint;

class CountHashGrid
{
public:
	~CountHashGrid()
	{
		clear();
	}

	void clear()
	{
		effectiveIndex.clear();
		effectiveIndex.shrink_to_fit();
		effectiveWeights.clear();
		effectiveWeights.shrink_to_fit();
	}

	void init(Scene *scene , int objIndex)
	{
		vec3f diag = scene->objKDTrees[objIndex].root->boundingBox.maxCoord -
			scene->objKDTrees[objIndex].root->boundingBox.minCoord;
		mBBoxMin = scene->objKDTrees[objIndex].root->boundingBox.minCoord;
		mBBoxMax = scene->objKDTrees[objIndex].root->boundingBox.maxCoord;

		sizeNum = 50;

		mCellSize = diag / sizeNum;
		mInvCellSize = 1.f / mCellSize;
		cellVolume = mCellSize.x * mCellSize.y * mCellSize.z;
		printf("cell size = (%.8f, %.8f, %.8f)\n" , mCellSize.x , mCellSize.y , mCellSize.z);
		
		sumWeights = 0.f;
		totVolume = 0.f;

		this->objectIndex = objIndex;
	}

	void preprocess(Scene *scene , int objIndex)
	{
		printf("preprocess volume object \#%d...\n" , objIndex);

		effectiveIndex.clear();
		effectiveWeights.clear();

		omp_lock_t lock;
		omp_init_lock(&lock);

#pragma omp parallel for
		for (int i = 0; i < sizeNum * sizeNum * sizeNum; i++)
		{
			Ray ray;
			float insideCnt = 0.f;
			int N = 8;
			for (int j = 0; j < N; j++)
			{
				vec3f offset;
				offset.x = RandGenerator::genFloat();
				offset.y = RandGenerator::genFloat();
				offset.z = RandGenerator::genFloat();
				ray.origin = cellIndexToPosition(i , offset);
				ray.direction = RandGenerator::genSphericalDirection();
				//ray.direction = vec3f(0.f , 1.f , 0.f);

				if (scene->checkInsideObject(ray , objIndex))
				{
					insideCnt += 1.f;
				}
			}
			float insideCellVol = insideCnt / (float)N;

			if (insideCellVol > 1e-6f)
			{
				omp_set_lock(&lock);
				sumWeights += insideCellVol;
				effectiveIndex.push_back(i);
				effectiveWeights.push_back(sumWeights);
				omp_unset_lock(&lock);
			}
		}

		for (int i = 0; i < effectiveWeights.size(); i++)
			effectiveWeights[i] /= sumWeights;

		totVolume = sumWeights * mCellSize.x * mCellSize.y * mCellSize.z;

		omp_destroy_lock(&lock);
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
			std::floor(mInvCellSize.x * distMin.x),
			std::floor(mInvCellSize.y * distMin.y),
			std::floor(mInvCellSize.z * distMin.z));

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

		vec3f corner = mBBoxMin + mCellSize * vec3f(x , y , z);

		vec3f res = corner + mCellSize * offset;

		return res;
	}

	vec3f getRandomPosition(float &pdf , int *givenIndex = NULL)
	{
		float rnd = RandGenerator::genFloat();
		unsigned index = (lower_bound(effectiveWeights.begin(), effectiveWeights.end(), rnd)-effectiveWeights.begin());

		if (givenIndex)
			index = *givenIndex;

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

	Ray emitVolume(Scene *scene , int *givenIndex = NULL)
	{
		for (;;)
		{
		Ray ray;
		ray.contactObject = NULL;
		ray.contactObjectTriangleID = -1;

		if (!givenIndex)
			ray.origin = getRandomPosition(ray.originProb);
		else
			ray.origin = getRandomPosition(ray.originProb , givenIndex);

		// not sure
		//ray.originProb /= cellVolume;
		ray.originProb = 1.f / totVolume;

		//printf("%.8f , %.8f\n" , ray.originProb , 1.f / scene->getTotalVolume());

		ray.direction = RandGenerator::genSphericalDirection();
		
		if (!scene->checkInsideObject(ray , objectIndex))
		{
			//printf("weird! It should always have insideObject\n");
			continue;
		}
		ray.insideObject = scene->objects[objectIndex];

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
			if(dist > 1e-6f)
			{
				ray.intersectDist = dist;
				ray.intersectObject = scene->objects[osi.objID];
				ray.intersectObjectTriangleID = osi.triangleID;
			}

			if (ray.insideObject && !ray.intersectObject)
			{
				printf("%.8f\n" , dist);
				printf("weird! It has insideObject but no intersectObject\n");
				continue;
			}
		}

		return ray;
		}
	}

	void print(FILE *fp , Scene* scene)
	{
		fprintf(fp , "============ one iter============\n");
		
		for (int i = 0; i < effectiveWeights.size(); i++)
		{
			Ray ray = emitVolume(scene , &i);
			fprintf(fp , "=====================\n");
			fprintf(fp , "pos = (%.8f, %.8f, %.8f) , dir = (%.8f, %.8f, %.8f)\n" , ray.origin.x , ray.origin.y , ray.origin.z ,
				ray.direction.x , ray.direction.y , ray.direction.z);
			fprintf(fp , "insideObj = %d, contactObj = %d, intersectObj = %d , t = %.8f\n" , ray.insideObject , ray.contactObject , 
				ray.intersectObject , ray.intersectDist);
			fprintf(fp , "index = %d, accuWeight = %.8f\n" , effectiveIndex[i] , effectiveWeights[i]);
		}
	}

public:
	vec3f mBBoxMin;
	vec3f mBBoxMax;

	std::vector<double> effectiveWeights;
	std::vector<int> effectiveIndex;

	vec3f mCellSize;
	vec3f mInvCellSize;

	float sumWeights;
	float cellVolume;
	float totVolume;
	int sizeNum;
	int objectIndex;
};

