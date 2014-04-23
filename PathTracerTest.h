#pragma once
#include "MCRenderer.h"

class PathTracerTest : public MCRenderer
{
private:
	unsigned spp;

	bool useConnection;

public:
	PathTracerTest(Renderer* renderer) : MCRenderer(renderer)
	{
		maxDepth = 20;
		spp = -1;
		useConnection = true;
	}

	bool mustUsePT(const Path& connectedPath)
	{
		if(connectedPath.back().contactObject && connectedPath.back().contactObject->emissive() &&
			(connectedPath.end()-2)->directionSampleType != Ray::RANDOM)
			return true;
		return false;
	}

	vec3f colorByConnectingLights(const Camera& camera, const Ray& ray, const Ray& lastRay);

	vector<vec3f> renderPixels(const Camera& camera);
};

