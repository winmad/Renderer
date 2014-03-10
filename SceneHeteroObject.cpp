#include "StdAfx.h"
#include "SceneHeteroObject.h"
#include <fstream>
#include <assert.h>


void SceneHeteroObject::parseDensityFile(const std::string &filename){
	std::ifstream fin(filename.c_str());
	fin >> nx >> ny >> nz;
	stepSize = 0.05;
	offSet = 0.00001;
	// TODO: parse the density file according to pbrt.
}

float SceneHeteroObject::getDensity(const vec3f &pos) const{
	// TODO: add code handling parsed density file, interpolation.
	return 1.0;
}

vec3f SceneHeteroObject::sigma_a(const vec3f &pos) const{
	return getDensity(pos) * (da);
}

vec3f SceneHeteroObject::sigma_s(const vec3f &pos) const{
	return getDensity(pos) * (ds);
}

vec3f SceneHeteroObject::sigma_t(const vec3f &pos) const{
	return getDensity(pos) * (dt);
}

float SceneHeteroObject::sampleDistance(const vec3f &pos, const vec3f &dir) const{
	float theta_t = y(dt);
	float t = -log(RandGenerator::genFloat()) / theta_t;
	while(y(sigma_t(pos + dir * t)) / theta_t < RandGenerator::genFloat()){
		t -= log(RandGenerator::genFloat()) / theta_t;
	}
	return t;
}

vec3f SceneHeteroObject::rayMarchingPdf(const vec3f &p0, const vec3f &p1, vec3f &averDt) const{
	vec3f direction = p1 - p0;
	float length = direction.length();
	if (length == 0.f) return 0.f;
	
	direction.normalize();
	float t = offSet;

	vec3f kernelPdf(0.,0.,0.);
	vec3f kernelDt(0.,0.,0.);
	int spp = 0;
	while(t < length){
		vec3f Dt = sigma_t(p0 + direction * t);
		kernelPdf += Dt * exp(-y(Dt) * t);
		kernelDt += Dt;
		t += stepSize;
		spp++;
	}
	assert(spp);
	averDt = kernelDt / spp;
	return kernelPdf * stepSize;
}

Ray SceneHeteroObject::scatter(Ray& inRay, float dist, unsigned triangleID, bool forward, SceneObject* intersectObject) const
{
	Ray outRay;
	if(forward)		return outRay;

	vec3f position = inRay.origin + inRay.direction*dist;
	LocalFrame lf = getAutoGenWorldLocalFrame(triangleID, position);

	if(intersectObject == this && inRay.direction.dot(lf.n) < 0)
	{
		outRay = inRay;
		outRay.origin = position;
		outRay.insideObject = (SceneObject*)this;
		outRay.contactObject = outRay.insideObject;
		return outRay;
	}
	outRay = inRay;
	float bounce_dist = sampleDistance(position, inRay.direction);
	bool hit = bounce_dist > dist;
	outRay.color = vec3f(1,1,1);
	outRay.directionProb = 1.0;
	outRay.insideObject = NULL;
	outRay.contactObject = NULL;

	if(hit){
		// hit a surface or volume bound 
		vec3f averDt;
		vec3f P_surf = P_surface(inRay.origin, inRay.origin + inRay.direction * dist, averDt);
		outRay.color = vec3f(1,1,1);
		outRay.origin = inRay.origin + inRay.direction * dist;
		outRay.direction = inRay.direction;
		outRay.directionProb *= y(P_surf);

		if(this == intersectObject)
		{
			outRay.insideObject = scene->findInsideObject(outRay, (SceneObject*)this);
			outRay.contactObject = NULL;
		}
		else
		{
			outRay.insideObject = (SceneObject*)this;
			outRay.contactObject = intersectObject;
		}
	}
	else{
		// simple travel in volume media 
		// travel straight forward
		outRay.color = vec3f(1, 1, 1);
		vec3f averDt;
		outRay.directionProb = y(p_medium(inRay.origin, inRay.origin + inRay.direction * bounce_dist, averDt));
		outRay.insideObject = (SceneObject*)this;
		outRay.contactObject = NULL;
		outRay.origin = inRay.origin + inRay.direction * bounce_dist;

		float albedo = y(ds) / y(dt);
		float rander = RandGenerator::genFloat();
		if(rander < albedo){
			// scattering 
			outRay.direction = vec3f(1,1,1);//RandGenerator::getHGDirection(-0.5, RandGenerator::genFloat(), RandGenerator::genFloat());
			vec3f phase_f = this->bsdf->evaluate(LocalFrame(), -outRay.direction, -inRay.direction);
			outRay.color *= phase_f;
			outRay.directionProb *= albedo * 1.0 / (4 * M_PI);
		}
		else{
			// terminate
			outRay.direction = vec3f(0, 0, 0); 
			outRay.color = vec3f(0, 0, 0);  
			outRay.directionProb = 1; 
			outRay.insideObject = NULL;
			outRay.contactObject = NULL;
		}
	}
 
	return outRay;
}
