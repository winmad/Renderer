#pragma once

#include <vector>
#include <cmath>
#include "nvVector.h"
using namespace nv;

typedef unsigned int uint;
typedef vec3f Vec3f;
typedef vec3i Vec3i;
typedef vec2i Vec2i;

class HashGrid
{
public:
	void Reserve(int aNumCells)
	{
		mCellEnds.resize(aNumCells);
	}

	template<typename tParticle>
	void Build(
		const std::vector<tParticle> &aParticles,
		float aRadius, float aMaxDist = 0.f)
	{
		if(std::fabsf(aMaxDist) < EPSILON)
			aMaxDist = aRadius;

		mMaxDistance = aMaxDist;
		mRadius      = aRadius;
		mCellSize    = aMaxDist * 2.f;
		mInvCellSize = 1.f / mCellSize;

		mBBoxMin = Vec3f( 1e36f);
		mBBoxMax = Vec3f(-1e36f);

		for(size_t i=0; i<aParticles.size(); i++)
		{
			const Vec3f &pos = aParticles[i].GetPosition();
			for(int j=0; j<3; j++)
			{
				mBBoxMax[j] = std::max(mBBoxMax[j], pos[j]);
				mBBoxMin[j] = std::min(mBBoxMin[j], pos[j]);
			}
		}

		mIndices.resize(aParticles.size());
		memset(&mCellEnds[0], 0, mCellEnds.size() * sizeof(int));

		// set mCellEnds[x] to number of particles within x
		for(size_t i=0; i<aParticles.size(); i++)
		{
			const Vec3f &pos = aParticles[i].GetPosition();
			mCellEnds[GetCellIndex(pos)]++;
		}

		// run exclusive prefix sum to really get the cell starts
		// mCellEnds[x] is now where the cell starts
		int sum = 0;
		for(size_t i=0; i<mCellEnds.size(); i++)
		{
			int temp = mCellEnds[i];
			mCellEnds[i] = sum;
			sum += temp;
		}

		for(size_t i=0; i<aParticles.size(); i++)
		{
			const Vec3f &pos = aParticles[i].GetPosition();
			const int targetIdx = mCellEnds[GetCellIndex(pos)]++;
			mIndices[targetIdx] = int(i);
		}

		// now mCellEnds[x] points to the index right after the last
		// element of cell x

		//// DEBUG
		//for(size_t i=0; i<aParticles.size(); i++)
		//{
		//    const Vec3f &pos  = aParticles[i].GetPosition();
		//    Vec2i range = GetCellRange(GetCellIndex(pos));
		//    bool found = false;
		//    for(;range.x < range.y; range.x++)
		//    {
		//        if(mIndices[range.x] == i)
		//            found = true;
		//    }
		//    if(!found)
		//        printf("Error at particle %d\n", i);
		//}
	}

	template<typename tParticle, typename tQuery>
	void Process(
		const std::vector<tParticle> &aParticles,
		tQuery& aQuery, bool useEllipse = false)
	{
		const Vec3f queryPos = aQuery.GetPosition();

		mMaxDistance = std::sqrt(mAniRatio) * mRadius;

		const Vec3f distMin = queryPos - mBBoxMin;
		const Vec3f distMax = mBBoxMax - queryPos;
		for(int i=0; i<3; i++)
		{
			if(distMin[i] < 0.f) return;
			if(distMax[i] < 0.f) return;
		}

		const Vec3f cellPt = mInvCellSize * distMin;
		const Vec3f coordF(
			std::floor(cellPt.x),
			std::floor(cellPt.y),
			std::floor(cellPt.z));

		const int  px = int(coordF.x);
		const int  py = int(coordF.y);
		const int  pz = int(coordF.z);

		const Vec3f fractCoord = cellPt - coordF;

		const int  pxo = px + (fractCoord.x < 0.5f ? -1 : +1);
		const int  pyo = py + (fractCoord.y < 0.5f ? -1 : +1);
		const int  pzo = pz + (fractCoord.z < 0.5f ? -1 : +1);

		//int foundCount = 0;

		for(int j=0; j<8; j++)
		{
			Vec2i activeRange;
			switch(j)
			{
			case 0: activeRange = GetCellRange(GetCellIndex(Vec3i(px , py , pz ))); break;
			case 1: activeRange = GetCellRange(GetCellIndex(Vec3i(px , py , pzo))); break;
			case 2: activeRange = GetCellRange(GetCellIndex(Vec3i(px , pyo, pz ))); break;
			case 3: activeRange = GetCellRange(GetCellIndex(Vec3i(px , pyo, pzo))); break;
			case 4: activeRange = GetCellRange(GetCellIndex(Vec3i(pxo, py , pz ))); break;
			case 5: activeRange = GetCellRange(GetCellIndex(Vec3i(pxo, py , pzo))); break;
			case 6: activeRange = GetCellRange(GetCellIndex(Vec3i(pxo, pyo, pz ))); break;
			case 7: activeRange = GetCellRange(GetCellIndex(Vec3i(pxo, pyo, pzo))); break;
			}

			for(; activeRange.x < activeRange.y; activeRange.x++)
			{
				const int particleIndex   = mIndices[activeRange.x];
				const tParticle &particle = aParticles[particleIndex];

				const float dist = (aQuery.GetPosition() - particle.GetPosition()).length();
				const float distSqr = dist * dist;


				if(distSqr <= mRadius * mRadius){
					aQuery.Process(particle);
				}

			}
		}
	}

private:

	Vec2i GetCellRange(int aCellIndex) const
	{
		if(aCellIndex == 0) return Vec2i(0, mCellEnds[0]);
		return Vec2i(mCellEnds[aCellIndex-1], mCellEnds[aCellIndex]);
	}

	int GetCellIndex(const Vec3i &aCoord) const
	{
		uint x = uint(aCoord.x);
		uint y = uint(aCoord.y);
		uint z = uint(aCoord.z);

		return int(((x * 73856093) ^ (y * 19349663) ^
			(z * 83492791)) % uint(mCellEnds.size()));
	}

	int GetCellIndex(const Vec3f &aPoint) const
	{
		const Vec3f distMin = aPoint - mBBoxMin;

		const Vec3f coordF(
			std::floor(mInvCellSize * distMin.x),
			std::floor(mInvCellSize * distMin.y),
			std::floor(mInvCellSize * distMin.z));

		const Vec3i coordI  = Vec3i(int(coordF.x), int(coordF.y), int(coordF.z));

		return GetCellIndex(coordI);
	}

private:
	Vec3f mBBoxMin;
	Vec3f mBBoxMax;
	std::vector<int> mIndices;
	std::vector<int> mCellEnds;

	float mMaxDistance;
	float mRadius;
	float mRadiusSqr;
	float mCellSize;
	float mInvCellSize;

public:
	vec3f		 mHitPoint;
	float		 mAniRatio;

};


