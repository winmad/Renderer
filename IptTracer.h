#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "PointKDTree.h"
#include "CountHashGrid.h"
#include "macros.h"

struct IptPathState
{
	vec3f throughput , indirContrib;
	//vec3f contribs[4]; 
	Ray *ray , *lastRay , *originRay;
	bool isSpecularPath;
	vec3f pos;
	int index;
	//double accuProb;
};

class IptTracer : public MCRenderer
{
protected:
	unsigned spp;

	int pixelNum , lightPathNum , cameraPathNum , interPathNum , partialPathNum;
	int lightPhotonNum , partialPhotonNum;
	int mergeIterations;

	vector<vector<int> > partPathMergeIndex;

	vector<float> weights;

	int *revIndex;

	vector<IptPathState> partialSubPathList;

	CountHashGrid countHashGrid;

	void movePaths(omp_lock_t& cmdLock , vector<Path>& , vector<Path*>&);

	void genLightPaths(omp_lock_t& cmdLock , vector<Path*>&);

	void genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>&);

	void mergePartialPaths(omp_lock_t& cmdLock);

	Ray genIntermediateSamples(vector<IptPathState>& partialSubPathList , Scene& scene , int *index = NULL);

	void calcEyeProbRatios(Path& eyePath , vector<float>& ratios);

	void mergePartialPaths(vector<vec3f>& contribs , const IptPathState& lightState , const int mergeIters);

	vec3f colorByMergingPaths(const IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths , const int mergeIters , Real& weight);

	vec3f colorByConnectingLights(const Camera& camera, const IptPathState& cameraState);

	vec3f colorByConnectingCamera(const Camera& camera, const IptPathState& lightState , int& _x , int& _y);

public:
	Real mergeRadius;
	Real mergeKernel;
	Real alpha;
	Real totArea , totVol;
	Real initialProb;
	unsigned timeInterval , lastTime;
	bool useWeight , usePPM;

	IptTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		alpha = 2.f / 3.f;
		spp = -1; 
		initialProb = 1.f;
		timeInterval = lastTime = 3600;

		pixelNum = renderer->camera.height * renderer->camera.width;
		cameraPathNum = pixelNum;
		lightPathNum = pixelNum;
		interPathNum = pixelNum;
		partialPathNum = pixelNum;

		usePPM = false;
		if (usePPM)
		{
			mergeIterations = 0;
			useWeight = false;
		}
		else
		{
			mergeIterations = 5;
			useWeight = true;
		}
	}
	void setRadius(const Real& r) { mergeRadius = r; }
	void setInitProb(const Real& r) { initialProb = r; }
	virtual vector<vec3f> renderPixels(const Camera& camera);
	
	Real getOriginProb(CountHashGrid& hashGrid , vec3f& pos , const bool isVol)
	{
		CountQuery query(pos);
		hashGrid.count(query);
		return query.count * (hashGrid.mInvCellSize * hashGrid.mInvCellSize * hashGrid.mInvCellSize);
	}
	
	
	Real connectFactor(Real pdf)
	{
		return pdf;
	}

	Real mergeFactor(Real *volScale = NULL , Real *initProb = NULL , Real *dirProb = NULL)
	{
		Real s = 1.0;
		if (volScale)
			s = *volScale;
		if (initProb)
			s *= *initProb;
		if (dirProb)
			s *= *dirProb;
		if (*volScale > 1.f)
			s *= 0.5f;
		Real res = M_PI * mergeRadius * mergeRadius * partialPathNum * s;
		return res;
	}

	Real calcEyePathWeight(Path& eyePath , vector<float>& ratios , int t)
	{
		Real sum = 1.f , tmp = 1.f;
		if (t - 1 >= ratios.size())
			return 1.f;
		for (int i = t - 1; i >= 0; i--)
		{
			tmp *= ratios[i];
			if (eyePath[i].directionSampleType == Ray::RANDOM)
				sum += tmp;
		}
		return 1.f / sum;
	};
};

