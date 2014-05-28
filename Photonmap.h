#pragma once
#include "MCRenderer.h"
#include "Hashgrid.h"
#include "BSDF.h"
#include "SceneVPMObject.h"
#include "SceneHeterogeneousVolume.h"
#include <omp.h>

class PhotonMap : public MCRenderer
{
private:
	HashGrid   surfaceHashGrid;
	omp_lock_t surfaceHashGridLock;
	HashGrid   volumeHashGrid;
	omp_lock_t volumeHashGridLock;
	float	   mBaseRadius;
	float	   mRadius;
	float	   mAlpha;
	int		   mPhotonsNum;
	unsigned   spp;
	unsigned   timeInterval , recordTime;
	bool	   rayMarching;

	typedef struct HashGridVertex{
		Path  *pathThePointIn;
		uint   indexInThePath;
		vec3f  position;
		Ray::PhotonType photonType;
		HashGridVertex(){
			pathThePointIn = NULL;
			photonType = Ray::OUTVOL;
			indexInThePath = -1;
		}
		vec3f GetPosition() const{ return position; }
	} LightPoint;
public:
	PhotonMap(Renderer *renderer) : MCRenderer(renderer)
	{
		mAlpha = 2.0/3;
		mPhotonsNum = 2 * renderer->camera.width * renderer->camera.height;
		spp = -1;
		timeInterval = 3600;
		recordTime = 3600;
	}
	void setRadius(float r) { mBaseRadius = r; }

	void throughputByDensityEstimation(vec3f &color, Path &eyeMergePath, 
		std::vector<LightPoint> &surfaceVertices, std::vector<LightPoint> &volumeVertices);

	void sampleMergePath(Path &path, Ray &prevRay, uint depth) const;

	virtual vector<vec3f> renderPixels(const Camera& camera);

	omp_lock_t debugPrintLock;
};



