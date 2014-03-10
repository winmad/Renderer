#include "StdAfx.h"
#include "Renderer.h"

Renderer::Renderer()
{
	camera.setScene(&scene);
	configManager = NULL;
	mcRenderer = NULL;
}

void Renderer::preview()
{
	IplImage* image = cvCreateImage(cvSize(camera.width, camera.height), IPL_DEPTH_32F, 3);
	vector<vec3f> eyeRays = camera.generateRays();
#pragma omp parallel for
	for(int x=0; x<camera.width; x++)
	{
		class SmoothNormal : public KDTree::Condition
		{
		public:
			Scene* scene;
			virtual bool legal(const KDTree::Ray& ray, const KDTree::Triangle& tri, const float dist) const
			{
				SceneObject *intersectObject = scene->objects[((Scene::ObjSourceInformation*)tri.sourceInformation)->objID];
				unsigned fi = ((Scene::ObjSourceInformation*)tri.sourceInformation)->triangleID;
				bool in = ray.direction.dot(intersectObject->getWorldNormal(fi, ray.origin + ray.direction*dist))<0;
				return in;
			}
		} condition;
		for(unsigned y=0; y<camera.height; y++)
		{
			condition.scene = &scene;
			Ray ray;
			ray.direction = eyeRays[y*camera.width + x];
			ray.origin = camera.position;
			Scene::ObjSourceInformation osi;
			float dist = scene.intersect(ray, osi, &condition);
			vec3f normal = dist >= 0 ? scene.objects[osi.objID]->getWorldNormal(osi.triangleID, ray.origin + ray.direction*dist) : vec3f(0, 0, 0);
			((vec3f*)image->imageData)[y*camera.width + x] = vec3f(1, 1, 1) * abs(ray.direction.dot(normal));
		}
	}
	cvShowImage("Renderer", image);
	cvReleaseImage(&image);
}

void Renderer::render()
{
	scene.buildObjKDTrees();
	vector<vec3f> pixelColors = mcRenderer->renderPixels(camera);
}

void Renderer::waitForCommand()
{
	int key = cvWaitKey(10);
	loadConfig("Data/Config.xml");
	render();
	/*
	while(key != 'q')
	{
		switch(key)
		{
		case 'p':
			loadConfig("Data/Config.xml");
			preview();
			break;
		case 'r':
			render();
			break;
		}
		key = cvWaitKey(10);
	}
	*/
}

void Renderer::loadConfig(const string& configFilePath)
{
	scene.clear();
	if(mcRenderer)
		delete mcRenderer;
	if(configManager)
		delete configManager;
	mcRenderer = NULL;
	configManager = NULL;
	configManager = new ConfigManager(this);
	configManager->load(configFilePath);
}


Renderer::~Renderer()
{
	if(mcRenderer)
		delete mcRenderer;
	if(configManager)
		delete configManager;
	mcRenderer = NULL;
	configManager = NULL;
	cvDestroyAllWindows();
}
