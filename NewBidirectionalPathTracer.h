#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "macros.h"

struct BidirPathState
{
	vec3f origin , pos;
	vec3f dir;
	vec3f throughput;
	BSDF bsdf;
	float dVCM , dVC;

	int pathLength : 15;
	int specularVertexNum : 15;
	int specularPath : 1;
	int isFiniteLight : 1;
};

class NewBidirectionalPathTracer : public MCRenderer
{
protected:
	unsigned spp;

	bool usePT;

	void genLightSample(Path& lightPath , BidirPathState& lightState);

	vec3f colorByConnectingCamera(const Camera& camera, const BidirPathState& lightState , const Ray& ray , const Ray& lastRay , int& _x , int& _y);

	void genCameraSample(const Camera& camera , Path& cameraPath , BidirPathState& cameraState);

	vec3f colorByHittingLight();

public:
	int controlLength;
	int lightPathNum , cameraPathNum , pixelNum;
	int iterations;

	std::vector<BidirPathState> lightStates;
	std::vector<BidirPathState> cameraStates;

	std::vector<int> lightStateIndex;
	std::vector<int> cameraStateIndex;

	unsigned timeInterval , lastTime;

	NewBidirectionalPathTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		pixelNum = renderer->camera.height * renderer->camera.width;
		lightPathNum = cameraPathNum = pixelNum;
		spp = -1;
		usePT = false;
		lastTime = timeInterval = 3600;
	}
	virtual vector<vec3f> renderPixels(const Camera& camera);

	float mis(const float& w)
	{
		return w;
	}
};