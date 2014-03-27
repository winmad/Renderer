#pragma once
#include "SceneObject.h"
#include "HGPhaseFunc.h"

class SceneHeteroObject : public SceneObject
{
public:
	bool isVolumeric() { return true; }
private:
	BSDF* bsdf;
	vec3f ds, da, dt;
	int nx, ny, nz;


	float y(vec3f c) const {
        const float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
	float *density;
	void  parseDensityFile(const std::string &filename);
	float getDensity(const vec3f &pos) const;
	vec3f sigma_a(const vec3f &pos) const;
	vec3f sigma_s(const vec3f &pos) const;
	vec3f sigma_t(const vec3f &pos) const;



	Ray scatter(Ray& inRay, float dist, unsigned triangleID, bool forward, SceneObject* intersectObject) const;
	

	float sampleDistance(const vec3f &p, const vec3f &d) const;

	float stepSize;
	float offSet;
	vec3f rayMarchingPdf(const vec3f &p0, const vec3f &p1, vec3f &averDt) const;


	// 1. p_medium 
	vec3f p_medium(const vec3f &p0, const vec3f &p1, vec3f &averDt) const{
		vec3f pm = rayMarchingPdf(p0, p1, averDt);
		return y(pm);
	}
	// 2. P_surface
	vec3f P_surface(const vec3f &p0, const vec3f &p1, vec3f &averDt) const{
		vec3f Pd = 1 - rayMarchingPdf(p0, p1, averDt);
		return Pd;
	}
 
public:
	SceneHeteroObject(Scene *scene) : SceneObject(scene) {
		//bsdf = new HGPhaseFunc;
		ds = vec3f(2.55, 3.21, 3.77);
		da = vec3f(0.0011, 0.0024, 0.014);
		dt = ds + da;
		canMerge = true;
	}
	~SceneHeteroObject(void){
		delete bsdf;
	}
};

