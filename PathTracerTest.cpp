#include "StdAfx.h"
#include "PathTracerTest.h"
#include "SceneEmissiveObject.h"

FILE *fp = fopen("test_pt.txt" , "w");

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
						float dirPdfA = eyePath[i].contactObject->getEmissionWeight() / eyePath[i].contactObject->totalArea;
						float mis = 1.f;
						if (i > 1 && !lastSpecular)
						{
							float cosine = eyePath[i].getContactNormal().dot(-eyePath[i - 1].direction);
							float dist = (eyePath[i].origin - eyePath[i - 1].origin).length();
							float dirPdfW = dirPdfA * dist * dist / abs(cosine);
							mis = lastPdfW / (lastPdfW + dirPdfW);
							/*
							fprintf(fp , "==================\n");
							fprintf(fp , "thr=(%.6f,%.6f,%.6f), contrib=(%.6f,%.6f,%.6f), pdfA=%.6f, pdfW=%.6f, lastPdfW=%.6f, cosine=%.6f, mis=%.6f\n" , 
								throughput[0] , throughput[1] , throughput[2] , contrib[0] ,
								contrib[1] , contrib[2] , dirPdfA , dirPdfW , lastPdfW , cosine , mis);
							*/
						}
						color += throughput * contrib * mis;
						break;
					}

					if (eyePath[i].directionSampleType == Ray::RANDOM)
					{
						vec3f contrib = colorByConnectingLights(camera , eyePath[i] , eyePath[i - 1]);
						color += throughput * contrib;
					}
					
					if (i + 1 >= eyePath.size() || eyePath[i].origin == eyePath[i + 1].origin)
						break;

					lastSpecular = (eyePath[i].directionSampleType == Ray::DEFINITE);
					lastPdfW = eyePath[i].directionProb;

					throughput *= eyePath[i].color * eyePath[i].getCosineTerm() / 
						eyePath[i].directionProb;
					float dist = std::max((eyePath[i].origin - eyePath[i - 1].origin).length() , 1e-5f);
					vec3f decayFactor = eyePath[i - 1].getRadianceDecay(dist);
					throughput *= decayFactor;
				}

				pixelColors[p] += color / ((float)s + 1);
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

	float cosAtLight = min(max(0.f , abs(lightRay.getContactNormal().dot(lightRay.direction))) , 1.f);
	if (cosAtLight < 1e-6f)
		return vec3f(0.f);

	float cosToLight = 1;
	if (ray.contactObject)
	{
		cosToLight = min(max(0.f , abs(ray.getContactNormal().dot(-lightRay.direction))) , 1.f);
	}

	float bsdfToLightPdf = inRay.getDirectionSampleProbDensity(outRay);
	outRay.direction = -lastRay.direction;
	float lightOriginPdf = lightRay.originProb;
	float dirPdfW = lightOriginPdf * dist2 / cosAtLight;

	float mis = dirPdfW / (dirPdfW + bsdfToLightPdf);

	vec3f tmp = lightRay.color * cosAtLight * bsdfFactor * cosToLight
		/ (lightRay.originProb * dist2);

	vec3f res = tmp * mis;
	/*
	fprintf(fp , "==================\n");
	fprintf(fp , "contrib=(%.6f,%.6f,%.6f), dirPdfW=%.6f, bsdfPdfW=%.6f, cosine=%.6f, mis=%.6f\n" , 
		lightRay.color[0] , lightRay.color[1] , lightRay.color[2] , dirPdfW , bsdfToLightPdf , cosToLight , mis);
	*/
	return res;
}