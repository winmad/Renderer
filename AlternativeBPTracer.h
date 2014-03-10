#pragma once
#include <omp.h>
#include "MCRenderer.h"
#include <unordered_map>
#include "macros.h"

class AlternativeBPTracer : public MCRenderer
{
public:
protected:
	unsigned spp;

	Ray link(const Ray& src, const Ray& dst);

	vec4<float> connectColorProb(const Path& connectedPath, int connectIndex);

	float connectWeight(int pixelID, int lightRaysCaught, const Path& connectedPath, int connectIndex, float expTerm = 1);

	void colorByConnectingPaths(vector<unordered_map<unsigned, unsigned>>& nLightRaysCaught, vector<omp_lock_t> &pixelLocks, const Camera& camera, vector<vec3f>& colors, const Path& eyePath, const Path& lightPath);

	bool connectRays(Path& path, int connectIndex, int pixelID);

	void statLightRaysCaught(vector<unordered_map<unsigned, unsigned>>& nLightRaysCaught, vector<omp_lock_t> &pixelLocks, const Path& lightPath);
public:
	AlternativeBPTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		maxDepth = 20;
		spp = -1;
	}
	virtual vector<vec3f> renderPixels(const Camera& camera);
};

