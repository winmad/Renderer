#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "macros.h"

class BidirectionalPathTracer : public MCRenderer
{
protected:
	unsigned spp;

	bool usePT;

	bool mustUsePT(const Path& connectedPath);

	Ray link(const Path& connectedPath, int connectIndex, int i, int j);

	vec4f connectColorProb(const Path& connectedPath, int connectIndex);

	float connectWeight(const Path& connectedPath, int connectIndex, vector<double>& p_forward, vector<double>& p_backward, vector<double>&distList, float expTerm = 1);

	void colorByConnectingPaths(vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const Path& eyePath, const Path& lightPath, vector<unsigned>* visibilityList = NULL);
public:
	unsigned timeInterval , lastTime;

	BidirectionalPathTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		spp = -1;
		usePT = false;
		lastTime = timeInterval = 3600;
	}
	virtual vector<vec3f> renderPixels(const Camera& camera);
};

