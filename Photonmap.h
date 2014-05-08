#pragma once
#include "MCRenderer.h"
#include "Hashgrid.h"
#include "BSDF.h"
#include "SceneVPMObject.h"

class Photonmap : public MCRenderer
{
private:
	HashGrid mHashgrid;
	HashGrid mVolHashgrid;
	unsigned spp;
	unsigned photonsWant;
	unsigned mResolution;

	float mBaseRadius;
	float mRadiusAlpha;
 
	float coneFilterK;
	vector<Ray> mLightVertices;
	vector<Ray> mVolLightVertices;

	void showPhotons();

	void shootPhotons();
	void buildHashgrid(unsigned pathCount, vector<Ray> &mLightVertices, float radius);
	void convLightpath2Photonpath(Path &lightPath, Path &photonPath, Path &volPhotonPath);
public:
	Photonmap(Renderer* renderer) : MCRenderer(renderer)
	{
		mRadiusAlpha = 0.75;
 
		photonsWant = renderer->camera.height * renderer->camera.width;
		mResolution = renderer->camera.height * renderer->camera.width;

		std::cout << "default photonsWant = " << photonsWant << std::endl;
		spp = 1000;
		coneFilterK = 1.0;
	}
	vector<vec3f> renderPixels(const Camera& camera);
	void setPhotonsWant(int num){
		photonsWant = num;
		std::cout << "set photonsWant = " << photonsWant << std::endl;
	}
	void setRadius(float r){
		mBaseRadius = r;
		std::cout << "set mBaseRadius = " << mBaseRadius << std::endl;
	}

	class RangeQuery{
		Ray &outRay;
		vec3f pos; 
		SceneObject *objPtr;
		bool isVol;
	public:
		RangeQuery(SceneObject *obj, Ray &outRay, vec3f pos, bool vol):objPtr(obj), outRay(outRay), pos(pos), isVol(vol){}
		void Process(const Ray &photon, vec3f &pixelColor){
			if(photon.photonType == Ray::OUTVOL && !isVol){
				pixelColor = (photon.color / photon.directionProb) * photon.getBSDF(outRay) / objPtr->getContinueProbability(photon, outRay); 
			}
			else if(photon.photonType == Ray::INVOL && isVol){
				pixelColor = photon.color / photon.directionProb * photon.getBSDF(outRay)  / objPtr->getContinueProbability(photon, outRay);
			}
		}
		vec3f GetPosition(){
			return pos;
		}
	};
	
	vec3f volMerge(SceneObject *obj, Ray &keyDiffRay);

	void sampleMergePath(Path& path, Ray& prevRay, unsigned depth) const;
};



