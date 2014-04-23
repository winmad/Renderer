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

	vec4f connectColorProb(const Path& connectedPath, int connectIndex);

	vector<vec3f> renderPixels(const Camera& camera);
};

