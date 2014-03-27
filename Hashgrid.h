#pragma once

#include <vector>
#include <cmath>
#include "nvVector.h"
using namespace nv;


typedef unsigned int uint;


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
        float aRadius, int aPhotonsWant)
    {
        mRadius      = aRadius;
        mRadiusSqr   = mRadius * mRadius;
        mCellSize    = mRadius * 2.f;
		mPhotonsWant = aPhotonsWant;
        mInvCellSize = 1.f / mCellSize;
        mBBoxMin = vec3f( 1e36f);
        mBBoxMax = vec3f(-1e36f);

        for(size_t i=0; i<aParticles.size(); i++)
        {
            const vec3f &pos = aParticles[i].origin;
            for(int j=0; j<3; j++)
            {
                mBBoxMax = std::max(mBBoxMax[j], pos[j]);
                mBBoxMin = std::min(mBBoxMin[j], pos[j]);
            }
        }

        mIndices.resize(aParticles.size());
        memset(&mCellEnds[0], 0, mCellEnds.size() * sizeof(int));

        // set mCellEnds[x] to number of particles within x
        for(size_t i=0; i<aParticles.size(); i++)
        {
            const vec3f &pos = aParticles[i].origin;
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
            const vec3f &pos = aParticles[i].origin;
            const int targetIdx = mCellEnds[GetCellIndex(pos)]++;
            mIndices[targetIdx] = int(i);
        }
 
    }

    template<typename tParticle, typename tQuery>
    void Process(
        const std::vector<tParticle> &aParticles,
        tQuery& aQuery, vec3f &color, bool isVolume)
    {
        const vec3f queryPos = aQuery.GetPosition();
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

        const int  px = int(coordF.x);
        const int  py = int(coordF.y);
        const int  pz = int(coordF.z);

        const vec3f fractCoord = cellPt - coordF;

        const int  pxo = px + (fractCoord.x < 0.5f ? -1 : +1);
        const int  pyo = py + (fractCoord.y < 0.5f ? -1 : +1);
        const int  pzo = pz + (fractCoord.z < 0.5f ? -1 : +1);

        int found = 0;

        for(int j=0; j<8; j++)
        {
            vec2i activeRange;
            switch(j)
            {
            case 0: activeRange = GetCellRange(GetCellIndex(vec3i(px , py , pz ))); break;
            case 1: activeRange = GetCellRange(GetCellIndex(vec3i(px , py , pzo))); break;
            case 2: activeRange = GetCellRange(GetCellIndex(vec3i(px , pyo, pz ))); break;
            case 3: activeRange = GetCellRange(GetCellIndex(vec3i(px , pyo, pzo))); break;
            case 4: activeRange = GetCellRange(GetCellIndex(vec3i(pxo, py , pz ))); break;
            case 5: activeRange = GetCellRange(GetCellIndex(vec3i(pxo, py , pzo))); break;
            case 6: activeRange = GetCellRange(GetCellIndex(vec3i(pxo, pyo, pz ))); break;
            case 7: activeRange = GetCellRange(GetCellIndex(vec3i(pxo, pyo, pzo))); break;
            }

            for(; activeRange.x < activeRange.y; activeRange.x++)
            {
                const int particleIndex   = mIndices[activeRange.x];
                const tParticle &particle = aParticles[particleIndex];

                const float aDist =
                    (aQuery.GetPosition() - particle.origin).length();
				const float distSqr = aDist * aDist;
                if(distSqr <= mRadiusSqr){
					vec3f colorCache = vec3f(0,0,0);
                    aQuery.Process(particle, colorCache);
					float kernel = Kernel(distSqr, mRadiusSqr);
					float weight = isVolume == false ? kernel / (mPhotonsWant * mRadiusSqr) : kernel / (mPhotonsWant * mRadiusSqr * mRadius);
					//float weight = isVolume == false ? 1 / (mPhotonsWant * M_PI * mRadiusSqr) : 1 / (mPhotonsWant * 4/3 * M_PI * mRadiusSqr * mRadius);
					color += colorCache *  weight;
				}
            }
        }
	}

private:

    vec2i GetCellRange(int aCellIndex) const
    {
        if(aCellIndex == 0) return vec2i(0, mCellEnds[0]);
        return vec2i(mCellEnds[aCellIndex-1], mCellEnds[aCellIndex]);
    }

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

private:

    vec3f mBBoxMin;
    vec3f mBBoxMax;
    std::vector<int> mIndices;
    std::vector<int> mCellEnds;
	int mPhotonsWant;

    float mRadius;
    float mRadiusSqr;
    float mCellSize;
    float mInvCellSize;
	float mConeFilter;
};

