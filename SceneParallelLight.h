#pragma once
#include "SceneEmissiveObject.h"

class SceneParallelLight : public SceneEmissiveObject
{
public:
	SceneParallelLight(Scene* scene) : SceneEmissiveObject(scene){
		canMerge = true;
	}
	Ray scatter(const Ray& inRay) const;
	Ray emit();
	float getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay);
};

