#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "PointKDTree.h"
#include "macros.h"
#include "smallFuncs.h"

struct LightPathPoint
{
	vec3f pos;
	Path* path;
	int index;
};


class VCMTracer : public MCRenderer
{
protected:
	unsigned spp;

	bool usePT;

	bool mustUsePT(const Path& connectedPath);

	Ray link(const Path& connectedPath, int connectIndex, int i, int j);

	vec4<float> connectColorProb(const Path& connectedPath, int connectIndex, bool merged = false);

	float connectMergeWeight(const Path& connectedPath, int connectIndex, bool merged, float expTerm = 1);

	void colorByMergingPaths(vector<vec3f>& colors, const Path& eyePath, PointKDTree<LightPathPoint>& tree);

	void colorByConnectingPaths(vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const Path& eyePath, const Path& lightPath, vector<unsigned>* visibilityList = NULL);

	
public:
	float mergeRadius;
	float alpha;
	unsigned timeInterval , lastTime;
	VCMTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		alpha = 0.6667;
		spp = -1;
		usePT = false;
		timeInterval = lastTime = 3600;
	}
	void setRadius(const float& r) { mergeRadius = r; }
	virtual vector<vec3f> renderPixels(const Camera& camera);
};

