#include "StdAfx.h"
#include "PathTracerTest.h"
#include "SceneEmissiveObject.h"

vector<vec3f> PathTracerTest::renderPixels(const Camera& camera)
{
	int t_start = clock();
	vector<vec3f> pixelColors(camera.width*camera.height, vec3f(0, 0, 0));

	if(useConnection)
		renderer->scene.preprocessEmissionSampler();

	if(!renderer->scene.usingGPU())
	{
		for(int s=0; s<spp; s++)
		{
			int t = clock();
#pragma omp parallel for
			for(int p=0; p<pixelColors.size(); p++)
			{
				Path eyePath;
				samplePath(eyePath, camera.generateRay(p));

				pixelColors[p] *= s/float(s+1);

				vec3f throughput = vec3f(1.f);
				vec3f color = vec3f(0.f);
				bool lastSpecular = 1;
				float lastPdfW = 1.f;

				for(unsigned i=1; i<eyePath.size(); i++)
				{
					if (eyePath[i].contactObject && eyePath[i].contactObject->emissive())
					{
						vec3f contrib = ((SceneEmissiveObject*)(eyePath[i].contactObject))->getColor();
						float dirPdfA = 1.f / eyePath[i].contactObject->getEmissionWeight();
						float mis = 1.f;
						if (i > 1 && !lastSpecular)
						{
							float cosine = eyePath[i].getContactNormal().dot(-eyePath[i - 1].direction);
							float dist = (eyePath[i].origin - eyePath[i - 1].origin).length();
							float dirPdfW = dirPdfA * dist * dist / cosine;
							mis = lastPdfW / (lastPdfW + dirPdfW);
						}
						color += throughput * contrib * mis;
						break;
					}

					vec3f contrib = colorByConnectingLights(camera , eyePath[i] , eyePath[i - 1]);
					color += throughput * contrib;

					if (i + 1 >= eyePath.size() || eyePath[i].origin == eyePath[i + 1].origin)
						break;

					lastSpecular = (eyePath[i].directionSampleType == Ray::DEFINITE);
					lastPdfW = eyePath[i].directionProb;

					throughput *= eyePath[i].color * eyePath[i].getCosineTerm() / 
						eyePath[i].directionProb;
				}

				pixelColors[p] += color;
			}

			showCurrentResult(pixelColors);
			printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, (clock()-t_start)/1000);
		}
	}
	return pixelColors;
}

vec3f PathTracerTest::colorByConnectingLights(const Camera& camera, const Ray& ray, const Ray& lastRay)
{
	Ray lightRay = genEmissiveSurfaceSample();
	lightRay.direction = (ray.origin - lightRay.origin);
	float dist = lightRay.direction.length();
	float dist2 = dist * dist;
	lightRay.direction.normalize();
	Ray outRay = ray;
	outRay.direction = -lightRay.direction;
	Ray inRay = lastRay;

	if(!testVisibility(outRay, lightRay))
		return vec3f(0.f);

	//outRay.direction = -cameraState.lastRay->direction;
	//vec3f bsdfFactor2 = lightRay.getBSDF(outRay);
	vec3f bsdfFactor = lastRay.getBSDF(outRay);

	if (y(bsdfFactor) < 1e-7f)
		return vec3f(0.f);

	float cosAtLight = min(max(0.f , lightRay.getContactNormal().dot(lightRay.direction)) , 1.f);

	float cosToLight = 1;
	if (ray.contactObject)
	{
		cosToLight = min(max(0.f , ray.getContactNormal().dot(-lightRay.direction)) , 1.f);
	}

	vec3f tmp = lightRay.color * cosAtLight * bsdfFactor * cosToLight
		/ (lightRay.originProb * dist2);

	float bsdfToLightPdf = inRay.getDirectionSampleProbDensity(outRay);

	outRay.direction = -lastRay.direction;
	float lightOriginPdf = lightRay.getOriginSampleProbDensity(outRay);
	float pdf = lightRay.getDirectionSampleProbDensity(outRay);

	float mis = lightOriginPdf * pdf / (lightOriginPdf * pdf + bsdfToLightPdf);

	vec3f res = tmp * mis;
	
	/*
	vec3f resx = camera.eliminateVignetting(res , cameraState.index) * pixelNum;
	if (cameraState.ray->contactObject)//(resx[0] + resx[1] + resx[2] >= 2)
	{
		fprintf(fp , "=====================\n");
		fprintf(fp , "cosAtLight = %.8f, cosToLight = %.8f, originPdf = %.8f, pdf = %.8f, weight=%.8f,\nres=(%.10f,%.10f,%.10f)\nbsdf=(%.10f,%.10f,%.10f)\n" , 
			cosAtLight , cosToLight , originProb , pdf , weightFactor , resx[0] , resx[1] , resx[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2]);
	}
	*/
	return res;
}