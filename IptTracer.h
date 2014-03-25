#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "PointKDTree.h"
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

	vector<vector<int> > partPathMergeIndex;

	vector<float> weights;

	int *revIndex;

	vector<IptPathState> partialSubPathList;

	void genLightPaths(omp_lock_t& cmdLock , vector<Path*>&);

	void genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>&);

	void mergePartialPaths(omp_lock_t& cmdLock);

	Ray genIntermediateSamples(vector<IptPathState>& partialSubPathList , Scene& scene);

	void mergePartialPaths(vector<vec3f>& contribs , const IptPathState& lightState);

	vec3f colorByMergingPaths(vector<vec3f>& colors, const IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths);

	vec3f colorByConnectingLights(const Camera& camera, vector<vec3f>& colors, const IptPathState& cameraState);

	void colorByConnectingCamera(vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const IptPathState& lightState);

public:
	Real mergeRadius;
	Real mergeKernel;
	Real alpha;
	Real totArea;
	Real initialProb;
	IptTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		maxDepth = 20;
		alpha = 0.75f;
		spp = -1; 
		initialProb = 1.f;

		pixelNum = lightPathNum = cameraPathNum = interPathNum = partialPathNum = 
			renderer->camera.height * renderer->camera.width;
	}
	void setRadius(const Real& r) { mergeRadius = r; }
	void setInitProb(const Real& r) { initialProb = r; }
	virtual vector<vec3f> renderPixels(const Camera& camera);

	Real connectFactor(Real pdf)
	{
		//return 0.5;
		return pdf;
	}

	Real mergeFactor(Real *volScale = NULL , Real *initProb = NULL , Real *dirProb = NULL)
	{
		//return 0.5;
		Real s = 1.0;
		if (volScale)
			s = *volScale;
		if (initProb)
			s *= *initProb;
		if (dirProb)
			s *= *dirProb;
		Real res = mergeRadius * mergeRadius * partialPathNum * s;
		return res;
	}
};

