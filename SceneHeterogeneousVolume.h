#pragma once
#include "SceneObject.h"
#include "HGPhaseSampler.h"
#include "HGPhaseFunc.h"
#define CLAMP(v, min_v, max_v) ((v < min_v) ? (min_v) : ((v > max_v) ? (max_v) : (v)))


class HeterogeneousVolume : public SceneObject{
	vec3f scatteringCoeff;
	vec3f absorptionCoeff;
	vec3f extinctionCoeff;
	HGPhaseFunc *bsdf;
	float g;
	typedef enum{ SCATTERING, ABSORPTION, EXTINCTION } LOOK_UP_TYPE; 
	float IOR;
	float stepSize; // for ray marching
	float densityScale;
	float scatteringScale;
	float absorptionScale;
	struct VolGridBox{
		vec3f p0, p1;
		int nx, ny, nz;
	} mBBox;
	float *mDensityMap;				// Smoke
	vec3f *mSubSurfaceDensityMap_Scattering;   // SubSurface
	vec3f *mSubSurfaceDensityMap_Absorption;   // SubSurface
	vec3f *mSubSurfaceDensityMap_Extinction;   // SubSurface
	inline float D(int x, int y, int z) const{
		x = CLAMP(x, 0, mBBox.nx-1);
		y = CLAMP(y, 0, mBBox.ny-1);
		z = CLAMP(z, 0, mBBox.nz-1);
		return mDensityMap[z*mBBox.nx*mBBox.ny + y*mBBox.nx + x] * densityScale;
	}
	bool  isSubsurface; // if subsurface, the density file will be RGB vec3f, instead of one float value.
	bool  volumeLookUpAll;
	inline vec3f SD(int x, int y, int z, LOOK_UP_TYPE type) const{
		x = CLAMP(x, 0, mBBox.nx-1);
		y = CLAMP(y, 0, mBBox.ny-1);
		z = CLAMP(z, 0, mBBox.nz-1);

		int index = z*mBBox.nx*mBBox.ny + y*mBBox.nx + x;
		vec3f densitySub;
		if(volumeLookUpAll){
			switch(type){
			case SCATTERING:
				densitySub = mSubSurfaceDensityMap_Scattering[index];	break;
			case ABSORPTION:
				densitySub = mSubSurfaceDensityMap_Absorption[index];	break;
			case EXTINCTION:
				densitySub = mSubSurfaceDensityMap_Extinction[index];	break;
			default:
				break;
			}
		}
		else{
			densitySub = mSubSurfaceDensityMap_Extinction[index];
		}
		return densitySub;
	}
	inline float lookUpDensity(const vec3f &worldPos) const;		   // Smoke
	inline vec3f lookUpSubSurfaceVolumeData(const vec3f &worldPos, LOOK_UP_TYPE type) const; // SubSurface
	inline bool checkIn(const vec3f &worldPos, const int objId) const;
	inline int check(const Ray &inRay, float *intersectDist = NULL) const;
public:
	HeterogeneousVolume(Scene *scene, float g) : SceneObject(scene),
		mDensityMap(NULL), 
		mSubSurfaceDensityMap_Scattering(NULL),
		mSubSurfaceDensityMap_Absorption(NULL),
		mSubSurfaceDensityMap_Extinction(NULL)
	{
		this->g = g;
		bsdf = new HGPhaseFunc(g);
		IOR = 1.5;
		stepSize = 0.1;
		densityScale = 100;
		scatteringScale = 100;
		absorptionScale = 100;
		canMerge = true;
		isSubsurface = false;
		volumeLookUpAll = true;
	}
	~HeterogeneousVolume(){
		delete bsdf;
		if(mDensityMap)					delete mDensityMap;
		if(mSubSurfaceDensityMap_Scattering)		delete mSubSurfaceDensityMap_Scattering;
		if(mSubSurfaceDensityMap_Absorption)		delete mSubSurfaceDensityMap_Absorption;
		if(mSubSurfaceDensityMap_Extinction)		delete mSubSurfaceDensityMap_Extinction;
	}
	void loadDensityMap(const std::string &filename);
	void loadSubSurfaceVolumeData(const std::string &filename, const std::string &filename2);
	void writeMitsubaDensityMap(const std::string& filename);
	void writeMitsubaAlbedo(const std::string& filename);

	Ray scatter(const Ray& inRay, const bool russian = true) const;
	vec3f getBSDF(const Ray &inRay, const Ray &outRay) const;
	float getOriginSampleProbDensity(const Ray &inRay, const Ray &outRay) const;
	float getDirectionSampleProbDensity(const Ray &inRay, const Ray &outRay) const;
	vec3f getRadianceDecay(const Ray &inRay, const float &dist) const;
	float getRefrCoeff() const { return IOR; }
	bool isVolumetric() { return true; }
	bool hasCosineTerm() { return false; }
	bool isHomogeneous() { return false; }
private:
	float integrateDensity(const Ray& inRay, float dist) const;
	vec3f tau(const Ray &inRay, float dist, bool noCheck = false) const;
	float pMedium(const Ray &inRay, float dist) const;
	float PSurface(const Ray &inRay, float dist) const;
	float getAlbedo() const;
	float getAlbedo(const vec3f &p) const;
	bool sampleDistance(const Ray &inRay, float &distance, float &pdfSuccess, float &pdfFailure) const;

	int findDesiredIntegralDensity(const Ray &inRay, const float desiredDensity, 
		float &t, float &integratedDensity, float &densityAtMinT, float &densityAtT) const;
public:
	void setSigma(vec3f sigmaS, vec3f sigmaA){
		scatteringCoeff = sigmaS;
		absorptionCoeff = sigmaA;
		extinctionCoeff = scatteringCoeff + absorptionCoeff;
		std::cout << "simga_s=" << scatteringCoeff << " sigma_a=" << absorptionCoeff << " simga_t= " << extinctionCoeff << std::endl;
	}
	void setDensityScale(float s){
		densityScale = s;
	}
	float getDensityScale() const{
		return densityScale;
	}
	void setScatteringScale(float s){
		scatteringScale = s;
	}
	float getScatteringScale() const{
		return scatteringScale;
	}
	void setAbsorptionScale(float s){
		absorptionScale = s;
	}
	float getAbsorptionScale() const{
		return absorptionScale;
	}
	void setStepSize(float s){
		stepSize = s;
	}
	float getStepSize() const{
		return stepSize;
	}
	void setIOR(float ior){
		IOR = ior;
		std::cout << "IOR=" << IOR << std::endl;
	}
	void setSubSurface(bool isSub){
		isSubsurface = isSub;
		return ;
	}
	void setVolumeLookUp(bool lookUpAll){
		// set 0: use all extinction coeffs.
		// set 1: use normal look up.
		volumeLookUpAll = lookUpAll;
		return ;
	}
};