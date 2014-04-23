#include "StdAfx.h"
#include "PathTracerTest.h"

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
						float dirPdfA = ((SceneEmissiveObject*)(eyePath[i].contactObject))->areaValues;
					}
				}

				pixelColors[p] += renderer->camera.eliminateVignetting(color, p)/(s+1)*camera.width*camera.height;
			}

			showCurrentResult(pixelColors);
			printf("Iter: %d  IterTime: %ds  TotalTime: %ds\n", s+1, (clock()-t)/1000, (clock()-t_start)/1000);
		}
	}
	return pixelColors;
}

