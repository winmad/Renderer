#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "PointKDTree.h"
#include "CountHashGrid.h"
#include "macros.h"

struct IptPathState
{
	vec3f throughput;
	vec3f dirContrib , indirContrib;
	Ray *ray , *lastRay , *originRay;
	bool isSpecularPath;
	vec3f pos;
	int index;
	double accuProb;
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

	void genLightPaths(omp_lock_t& cmdLock , vector<Path*>&);

	void genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>&);

	void mergePartialPaths(omp_lock_t& cmdLock);

	Ray genIntermediateSamples(vector<IptPathState>& partialSubPathList , Scene& scene);

	void mergePartialPaths(vector<vec3f>& contribs , const IptPathState& lightState);

	void calcEyeProbRatios(Path& eyePath , vector<float>& ratios);

	vec3f colorByMergingPaths(const IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths);

	vec3f colorByConnectingLights(const Camera& camera, const IptPathState& cameraState);

	vec3f colorByConnectingCamera(const Camera& camera, const IptPathState& lightState , int& _x , int& _y);

public:
	Real mergeRadius;
	Real mergeKernel;
	Real alpha;
	Real totArea , totVol;
	Real initialProb;
	unsigned timeInterval , lastTime;
	IptTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		maxDepth = 20;
		alpha = 0.75f;
		spp = -1; 
		initialProb = 1.f;
		mergeIterations = 5;
		timeInterval = lastTime = 3600;

		pixelNum = renderer->camera.height * renderer->camera.width;
		cameraPathNum = pixelNum;
		lightPathNum = pixelNum;
		interPathNum = pixelNum;
		partialPathNum = pixelNum;
	}
	void setRadius(const Real& r) { mergeRadius = r; }
	void setInitProb(const Real& r) { initialProb = r; }
	virtual vector<vec3f> renderPixels(const Camera& camera);

	Real getOriginProb(CountHashGrid& hashGrid , vec3f& pos , const bool isVol)
	{
		CountQuery query(pos);
		hashGrid.count(query);
		if (isVol)
			return query.count / hashGrid.cellVolume;
		else
			return query.count / hashGrid.cellArea;
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

