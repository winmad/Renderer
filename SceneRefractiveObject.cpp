#include "StdAfx.h"
#include "SceneRefractiveObject.h"

Ray SceneRefractiveObject::scatter(Ray& inRay) const
{

	Ray outRay;
	vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;

	LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, position);
	vec3f normal = lf.n;

	// scatter--start

	if(inRay.intersectObject != this && inRay.insideObject == this)
	{
		outRay = inRay.intersectObject->scatter(inRay);
		//outRay.directionSampleType = Ray::DEFINITE;
		return outRay;
	}

	outRay.origin = position;
	outRay.direction = inRay.direction;
	vec3f reflDir = -normal.dot(inRay.direction)*normal*2 + inRay.direction;
	reflDir.normalize();
	float theta = acos(inRay.direction.dot(normal));

	SceneObject* currentInsideObject = inRay.insideObject;
	SceneObject* outSideObject = (SceneObject*)this;
	if(inRay.insideObject == this)
		outSideObject = scene->findInsideObject(outRay, (SceneObject*)this);
	float current_n = currentInsideObject ? currentInsideObject->getRefrCoeff() : 1;
	float next_n = outSideObject ? outSideObject->getRefrCoeff() : 1;
	float sin_phi = current_n / next_n * sin(theta);

	outRay.intersectObject = NULL;
	outRay.color = surfColor;
	outRay.directionProb = 1;
	outRay.contactObject = (SceneObject*)this;
	outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;

	outRay.directionSampleType = Ray::DEFINITE;

	if(sin_phi > 1)
	{
		outRay.direction = reflDir;
		outRay.insideObject = inRay.insideObject;
		outRay.directionProb = 1;
		outRay.color /= outRay.getCosineTerm();
	}
	else
	{
		float phi = asin(sin_phi);
		if(theta > M_PI/2)
			phi = M_PI - phi;
		vec3f axis = normal.cross(inRay.direction);
		axis.normalize();
		outRay.direction = vec3f(rotMat(axis, phi) * vec4<float>(normal, 0));
		outRay.direction.normalize();

		float cos_theta = abs(cos(theta));
		float cos_phi = abs(cos(phi));
		float esr = powf(abs(current_n*cos_theta-next_n*cos_phi)/(current_n*cos_theta+next_n*cos_phi),2);
		float epr = powf(abs(next_n*cos_theta-current_n*cos_phi)/(next_n*cos_theta+current_n*cos_phi),2);
		float er = (esr+epr)/2;

		float p = er;

		if(RandGenerator::genFloat() < p)
		{
			outRay.direction = reflDir;
			outRay.color *= er / outRay.getCosineTerm();
			outRay.photonProb = outRay.directionProb = p;
			outRay.insideObject = inRay.insideObject;
		}
		else
		{
			outRay.color *= (1-er) / outRay.getCosineTerm();
			outRay.photonProb = outRay.directionProb = 1-p;
			outRay.insideObject = outSideObject;
		}
	}
	outRay.direction.normalize();
	
	return outRay;
}
