#pragma once
#include "SceneObject.h"
#include "nvVector.h"
#include <vector>

using namespace nv;

using namespace std;



class Camera : public SceneObject
{
public:
	vec3f focus;
	vec3f position;
	vec3f up;
	unsigned sightDist;
	int width;
	int height;
	void setScene(Scene* scene){ this->scene = scene; }
	float get_dw(unsigned pixelID) const;
	Ray generateRay(unsigned pixelID , bool flag = false) const;
	vec3f eliminateVignetting(const vec3f& color, unsigned pixelID) const;
	vector<vec3f> generateRays() const;
	vec2<float> transToPixel(const vec3f& point) const;
	virtual float getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const;
	virtual vec3f getBSDF(const Ray& inRay, const Ray& outRay) const { return vec3f(1, 1, 1); }
	vec3f getWorldNormal(unsigned fi, const vec3f& position, bool flat = false) const;
};

