#pragma once

#include <vector>
#include <cmath>
#include "nvVector.h"
using namespace nv;

typedef unsigned int uint;

struct CountQuery
{
	vec3f& pos;
	float count;

	CountQuery(vec3f& _pos) : pos(_pos) , count(0) {}
};

class CountHashGrid
{
public:
	void clear()
	{
		mCellEnds.clear();
		mCellEnds.shrink_to_fit();
	}

	void reserve(int aNumCells)
	{
		mCellEnds.resize(aNumCells);
	}

	template<typename tParticle>
	void build(const std::vector<tParticle> &aParticles,
		float aRadius)
	{
		mBBoxMin = vec3f( 1e36f);
		mBBoxMax = vec3f(-1e36f);

		for(size_t i=0; i<aParticles.size(); i++)
		{
			const vec3f &pos = aParticles[i].pos;
			for(int j=0; j<3; j++)
			{
				mBBoxMax = std::max(mBBoxMax[j], pos[j]);
				mBBoxMin = std::min(mBBoxMin[j], pos[j]);
			}
		}

		mCellSize = aRadius * 2.f;

		mCellEnds.resize(aParticles.size() / 5);
		memset(&mCellEnds[0], 0, mCellEnds.size() * sizeof(int));

		// set mCellEnds[x] to number of particles within x
		for(size_t i=0; i<aParticles.size(); i++)
		{
			const vec3f &pos = aParticles[i].pos;
			vec3f totContrib = aParticles[i].dirContrib + aParticles[i].indirContrib;
			mCellEnds[GetCellIndex(pos)] += y(totContrib);
		}

		sumContribs = 0;
		for (size_t i = 0; i < mCellEnds.size(); i++)
			sumContribs += mCellEnds[i];
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

		const int  pxo = px + (fractCoord.x < 0.5f ? -1 : +1);
		const int  pyo = py + (fractCoord.y < 0.5f ? -1 : +1);
		const int  pzo = pz + (fractCoord.z < 0.5f ? -1 : +1);

		int found = 0;

		for(int j=0; j<8; j++)
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

			aQuery.count += mCellEnds[cellIndex];
		}
	}

public:
	float Kernel(float distSqr, float radiusSqr) const{
		float s = 1 - distSqr / radiusSqr;
		return 3 * s * s / M_PI;
	}

	int GetCellIndex(const vec3i &aCoord) const
	{
		uint x = uint(aCoord.x);
		uint y = uint(aCoord.y);
		uint z = uint(aCoord.z);

		return int(((x * 73856093) ^ (y * 19349663) ^
			(z * 83492791)) % uint(mCellEnds.size()));
	}

	int GetCellIndex(const vec3f &aPoint) const
	{
		const vec3f distMin = aPoint - mBBoxMin;

		const vec3f coordF(
			std::floor(mInvCellSize * distMin.x),
			std::floor(mInvCellSize * distMin.y),
			std::floor(mInvCellSize * distMin.z));

		const vec3i coordI  = vec3i(int(coordF.x), int(coordF.y), int(coordF.z));

		return GetCellIndex(coordI);
	}

public:
	vec3f mBBoxMin;
	vec3f mBBoxMax;
	std::vector<float> mCellEnds;

	float mCellSize;
	float mInvCellSize;

	float sumContribs;
};

