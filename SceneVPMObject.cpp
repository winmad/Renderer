#include "StdAfx.h"
#include "SceneVPMObject.h"

Ray SceneVPMObject::scatter(const Ray& inRay, const bool russian) const
{
	Ray outRay;
	outRay.directionSampleType = Ray::RANDOM;
	HomoMediaDistSampler logDistSampler(y(dt));
	float sampDist = logDistSampler.sampleDist();

	bool go_in_vol  = (inRay.intersectObject == this) && (inRay.insideObject != this);
	bool be_in_vol  = (inRay.insideObject == this);
	bool out_of_vol = (sampDist >= inRay.intersectDist); // hit a surface

	// go in volume
	if(go_in_vol){
		vec3f position = inRay.origin + inRay.direction*inRay.intersectDist;
		LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, position);
		vec3f normal = lf.n;
		outRay.origin = position;
		outRay.direction = inRay.direction;
		
		vec3f reflDir = -normal.dot(inRay.direction)*normal*2 + inRay.direction;
		reflDir.normalize();
		float theta = acos(clampf(inRay.direction.dot(normal), -1, 1));
		
		SceneObject* currentInsideObject = inRay.insideObject;
		SceneObject* outSideObject = (SceneObject*)this;

		float current_n = currentInsideObject ? currentInsideObject->getRefrCoeff() : 1;
		float next_n = outSideObject ? outSideObject->getRefrCoeff() : 1;
		float sin_phi = current_n / next_n * sin(theta);

		outRay.intersectObject = NULL;
		outRay.color = vec3f(1, 1, 1);
		outRay.directionProb = 1;
		outRay.contactObject = (SceneObject*)this;
		outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;

		if(sin_phi > 1){
			outRay.direction = reflDir;
			outRay.insideObject = inRay.insideObject;
			outRay.directionProb = 1;
			outRay.directionSampleType = Ray::DEFINITE;
			outRay.photonType = Ray::NOUSE;
		}
		else{
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
				outRay.directionProb = p;
				outRay.insideObject = inRay.insideObject;
				outRay.directionSampleType = Ray::DEFINITE;
				outRay.photonType = Ray::NOUSE;
			}
			else
			{
				outRay.color *= (1-er) / outRay.getCosineTerm();
				outRay.directionProb = 1-p;
				outRay.contactObject = outRay.insideObject = (SceneObject*)this;
				outRay.directionSampleType = Ray::DEFINITE;
				outRay.photonType = Ray::HITVOL;
				outRay.isDirectLightPhoton = inRay.isDirectLightPhoton;
			}
			outRay.direction.normalize();
		}

		return outRay;
	}
	// in volume
	if(be_in_vol && !out_of_vol){
		HGPhaseSampler hgPhaseSampler(g);
		//IsotropicPhaseSampler sp;
		outRay.origin = inRay.origin + inRay.direction * sampDist;
		LocalFrame lf;		lf.n = inRay.direction;
		outRay.direction = hgPhaseSampler.genSample(lf); 
		outRay.insideObject = (SceneObject*)this;
		outRay.contactObject = NULL;
		//outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;
		outRay.contactObjectTriangleID = -1;

		outRay.directionProb = 1;
		outRay.color = vec3f(1, 1, 1);
		float albedo = y(ds) / y(dt);
		float rander = RandGenerator::genFloat();
		outRay.originSampleType = Ray::RANDOM;

		if(rander < albedo || (!russian)){
			outRay.contactObject = NULL;
			LocalFrame lf;	lf.n = inRay.direction;
			float oPdfW = hgPhaseSampler.getProbDensity(lf, outRay.direction);//sp.getProbDensity(lf, outRay.direction);//
			outRay.directionProb *= albedo * oPdfW;
			outRay.originProb = p_medium(sampDist);
			outRay.directionSampleType = Ray::RANDOM;
			outRay.color *= ds * bsdf->evaluate(LocalFrame(), inRay.direction, outRay.direction);
			outRay.photonType = inRay.isDirectLightPhoton ? Ray::NOUSE : Ray::INVOL;
			outRay.isDirectLightPhoton = false;
		}
		else{
			// terminate
			outRay.direction = vec3f(0, 0, 0); 
			outRay.color = vec3f(0, 0, 0);  
			outRay.directionProb = 1; 
			outRay.originProb = p_medium(sampDist);
			
			outRay.insideObject = (SceneObject*)this; // FIXED
			
			outRay.contactObject = NULL;
			outRay.directionSampleType = Ray::RANDOM;
			outRay.photonType = inRay.isDirectLightPhoton ? Ray::NOUSE : Ray::INVOL;
			outRay.isDirectLightPhoton = false;
		}
		return outRay;
	}
	// go out of volume
	if(be_in_vol && out_of_vol){
		outRay = inRay;
		outRay.direction = inRay.direction;
		outRay.origin = inRay.origin + inRay.intersectDist * inRay.direction;
		outRay.contactObject = inRay.intersectObject;
		outRay.contactObjectTriangleID = inRay.intersectObjectTriangleID;
		outRay.insideObject = (SceneObject*)this;
		outRay.directionProb = 1; 
		outRay.color = vec3f(1,1,1);
		bool going_out = (inRay.intersectObject == this);
		
		if(going_out){
			LocalFrame lf = inRay.intersectObject->getAutoGenWorldLocalFrame(inRay.intersectObjectTriangleID, outRay.origin);
			vec3f normal = lf.n;
			vec3f reflDir = -normal.dot(inRay.direction)*normal*2 + inRay.direction;
			reflDir.normalize();
			float theta = acos(clampf(inRay.direction.dot(normal), -1, 1));
			
			SceneObject* currentInsideObject = (SceneObject*)this;
			SceneObject* outSideObject = scene->findInsideObject(outRay, (SceneObject*)this);

			float current_n = currentInsideObject ? currentInsideObject->getRefrCoeff() : 1;
			float next_n = outSideObject ? outSideObject->getRefrCoeff() : 1;
			float sin_phi = current_n / next_n * sin(theta);

			outRay.intersectObject = NULL;
			if(sin_phi > 1){
				outRay.direction = reflDir;
				outRay.insideObject = inRay.insideObject;
				outRay.contactObject = (SceneObject*)this;
				outRay.originProb = P_surface(inRay.intersectDist);
				outRay.photonType = Ray::NOUSE;
				outRay.directionSampleType = Ray::DEFINITE;
			}
			else{
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
					outRay.directionProb = p;
					outRay.originProb =  P_surface(inRay.intersectDist);
					outRay.insideObject = inRay.insideObject;
					outRay.directionSampleType = Ray::DEFINITE;
					outRay.photonType = Ray::NOUSE;
				}
				else
				{
					outRay.color *= (1-er) / outRay.getCosineTerm();
					outRay.directionProb = (1-p);
					outRay.originProb = P_surface(inRay.intersectDist);
					outRay.insideObject = outSideObject;
					outRay.directionSampleType = Ray::DEFINITE;
					outRay.photonType = Ray::NOUSE;
				}
				outRay.direction.normalize();
			}	
		}
		else{
			//std::cout << "in vol hit surface " << std::endl;
			outRay.contactObject = NULL;
			outRay.intersectDist = 0;
			
			//------ FIX ME -------
			if (inRay.intersectObject != NULL)
				outRay = inRay.intersectObject->scatter(outRay , russian);
			else
				printf("error\n");
			//---------------------

			outRay.originProb *= P_surface(inRay.intersectDist);
			outRay.photonType = Ray::NOUSE;
		}
		return outRay;
	}
 
}

vec3f SceneVPMObject::getRadianceDecay(const Ray &inRay, const float &dist) const {
	return transmittance(dist);
}
 
vec3f SceneVPMObject::getBSDF(const Ray& inRay, const Ray& outRay) const
{
	LocalFrame lf;
	lf.n = inRay.direction;
	if(outRay.contactObject == NULL)
		return ds * bsdf->evaluate(LocalFrame(), inRay.direction, outRay.direction);
	if(outRay.contactObject && outRay.contactObject != this)
		return outRay.contactObject->getBSDF(inRay, outRay);
	return vec3f(0, 0, 0);
}

float SceneVPMObject::getOriginSampleProbDensity(const Ray &inRay, const Ray &outRay) const{
	
	float dist = max((inRay.origin - outRay.origin).length(), EPSILON);

	if(outRay.contactObject)
		return P_surface(dist);
	else
		return p_medium(dist);
	return 1;
}

float SceneVPMObject::getDirectionSampleProbDensity(const Ray &inRay, const Ray &outRay) const{
	if(outRay.directionSampleType == Ray::DEFINITE)
	{
		return 0;
	}
	if(outRay.contactObject)
		return outRay.contactObject->getDirectionSampleProbDensity(inRay, outRay);
	HGPhaseSampler hgPhaseSampler(g);
	LocalFrame lf;	lf.n = inRay.direction;
	float albedo = y(ds) / y(dt);
	float oPdfW = hgPhaseSampler.getProbDensity(lf, outRay.direction);
	return albedo*oPdfW;
}

float SceneVPMObject::getContinueProbability(const Ray &inRay, const Ray &outRay) const {
	float albedo = y(ds) / y(dt);
	return albedo;
}