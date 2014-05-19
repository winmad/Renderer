#include "StdAfx.h"
#include "MCRenderer.h"
#include "NoSelfIntersectionCondition.h"

void mouseEvent(int evt, int x, int y, int flags, void* param)
{
	MCRenderer *mcRenderer = (MCRenderer*)param;
	
	if(evt == CV_EVENT_MBUTTONDOWN)
	{ 
		mcRenderer->showEvent = MCRenderer::SHOW_PATH;
		mcRenderer->pathPixelID = y * mcRenderer->renderer->getCamera().width + x;
	}
	if(evt == CV_EVENT_LBUTTONDOWN)
	{
		mcRenderer->showEvent = MCRenderer::SHOW_PIXEL_VALUE;
		mcRenderer->pathPixelID = y * mcRenderer->renderer->getCamera().width + x;
	}
	if(evt == CV_EVENT_RBUTTONDOWN)
	{
		mcRenderer->showEvent = MCRenderer::NO_EVENT;
	}
}

void MCRenderer::setSavePath(const string& savePath)
{
	unsigned pos1 = savePath.rfind('/');
	unsigned pos2 = savePath.rfind('\\');
	unsigned pos = pos1 < pos2 ? pos1 : pos2;
	string folderPath = savePath.substr(0, pos);
	if(_access(folderPath.c_str(), 0))
	{
		this->savePath = renderer->configManager->getRootPath() + string("/") + savePath;
	}
}

void MCRenderer::resetInsideObject(Ray& ray)
{
	if(!ray.contactObject)
		return;
	vec3f normal = ray.contactObject->getWorldNormal(ray.contactObjectTriangleID, ray.origin);
	if(ray.direction.dot(normal)<0)
		ray.insideObject = ray.contactObject;
	else
		ray.insideObject = renderer->scene.findInsideObject(ray, ray.contactObject);
}

bool MCRenderer::connectRays(Path& path, int connectIndex, bool merged)
{
	Ray& lightRay = path[connectIndex];
	Ray& eyeRay = path[connectIndex+1];

	if(!merged && (lightRay.directionSampleType == Ray::DEFINITE || eyeRay.directionSampleType == Ray::DEFINITE))
		return false;

	if (merged && eyeRay.directionSampleType == Ray::DEFINITE)
		return false;

	const Ray& prevLightRay = path[max2(connectIndex - 1, 0)];
	const Ray& prevEyeRay = path[min2((connectIndex + 2), path.size()-1)];
	lightRay.direction = eyeRay.origin - lightRay.origin;
	lightRay.direction.normalize();
	eyeRay.direction = -lightRay.direction;

	if(!merged)
	{
		lightRay.directionProb = 1;
		lightRay.color = prevLightRay.getBSDF(lightRay);
		lightRay.directionSampleType = Ray::RANDOM;
		if(lightRay.contactObject && connectIndex > 0)
		{
			if(lightRay.getContactNormal().dot(lightRay.direction) <= 0)
			{
				lightRay.insideObject = lightRay.contactObject;
			}
			else
			{
				if(lightRay.getContactNormal().dot(prevLightRay.direction) <= 0)
					lightRay.insideObject = prevLightRay.insideObject;
				else
					lightRay.insideObject = renderer->scene.findInsideObject(lightRay, lightRay.contactObject);
			}
		}
	}

	eyeRay.directionProb = 1;
	eyeRay.color = prevEyeRay.getBSDF(eyeRay);
	eyeRay.directionSampleType = Ray::RANDOM;
	if(eyeRay.contactObject && connectIndex < path.size()-2)
	{
		if(eyeRay.getContactNormal().dot(eyeRay.direction) <= 0)
		{
			eyeRay.insideObject = eyeRay.contactObject;
		}
		else
		{
			if(eyeRay.getContactNormal().dot(prevEyeRay.direction) <= 0)
				eyeRay.insideObject = prevEyeRay.insideObject;
			else
				eyeRay.insideObject = renderer->scene.findInsideObject(eyeRay, eyeRay.contactObject);
		}
	}
	return true;
}

bool MCRenderer::testVisibility(const Ray& startRay, const Ray& endRay)
{
	/*NoSelfIntersectionCondition srCond(&renderer->scene, startRay);
	IntersectInfo srInfo = intersect(startRay.origin, startRay.direction, &srCond);
	NoSelfIntersectionCondition erCond(&renderer->scene, endRay);
	IntersectInfo erInfo = intersect(endRay.origin, endRay.direction, &erCond);
	if(srInfo.dist<0 || erInfo.dist<0)
		return true;
	return (srInfo.dist + erInfo.dist) > 1.5 * (startRay.origin - endRay.origin).length();*/
	NoSelfIntersectionCondition srCond(&renderer->scene, startRay);
	IntersectInfo srInfo = intersect(startRay.origin, startRay.direction, &srCond);
	return srInfo.dist<0 || (srInfo.dist - (startRay.origin - endRay.origin).length() >= -EPSILON);
}

SceneObject* MCRenderer::findInsideObject(const vec3f& origin, const vec3f& direction)
{
	Ray ray;
	ray.direction = direction;
	ray.origin = origin;
	return renderer->scene.findInsideObject(ray);
}

void MCRenderer::eliminateVignetting(vector<vec3f>& pixelColors)
{
	for(unsigned i=0; i<pixelColors.size(); i++)
	{
		pixelColors[i] = renderer->camera.eliminateVignetting(pixelColors[i], i);
	}
}

MCRenderer::IntersectInfo MCRenderer::intersect(const vec3f& origin, const vec3f& direction, KDTree::Condition* condition)
{
	Ray ray;
	ray.direction = direction;
	ray.origin = origin;
	IntersectInfo info;
	Scene::ObjSourceInformation osi;
	info.dist = renderer->scene.intersect(ray, osi, condition);
	info.triangleID = osi.triangleID;
	info.intersectObject = (info.dist >= 0) ? renderer->scene.objects[osi.objID] : NULL;
	return info;
}

void MCRenderer::samplePath(Path& path, Ray& prevRay, unsigned depth, bool firstDiff) const
{
	if(prevRay.directionSampleType == Ray::RANDOM && firstDiff){
		path.push_back(prevRay);
		return ;
	}
	Ray termRay;
	termRay.origin = prevRay.origin;
	termRay.direction = vec3f(0, 0, 0);
	termRay.color = vec3f(0, 0, 0);
	termRay.directionProb = 1;
	termRay.directionSampleType = Ray::DEFINITE;
	termRay.insideObject = termRay.contactObject = termRay.intersectObject = NULL;

	path.push_back(prevRay);

	Ray nextRay;

	if(prevRay.insideObject)
	{
		nextRay = prevRay.insideObject->scatter(prevRay);
	}
	else if(prevRay.intersectObject)
	{
		nextRay = prevRay.intersectObject->scatter(prevRay);
	}
	else
	{
		path.push_back(termRay);
		return;
	}

	if(nextRay.direction.length() < 0.5)
	{
		path.push_back(nextRay);
		return;
	}
	
	if(depth+1 > maxDepth)
	{
		path.push_back(termRay);
		return;
	}

	NoSelfIntersectionCondition condition(&renderer->scene, nextRay);

	Scene::ObjSourceInformation osi;
	float dist;
	dist = renderer->scene.intersect(nextRay, osi, &condition);

	if(dist < 0)
	{
		path.push_back(nextRay);
		path.push_back(termRay);
		return;
	}
	else
	{
		nextRay.intersectObject = renderer->scene.objects[osi.objID];
		nextRay.intersectObjectTriangleID = osi.triangleID;
		nextRay.intersectDist = dist;
	}
	
	samplePath(path, nextRay, depth + 1, firstDiff);
}

void MCRenderer::samplePath(Path& path, Ray& startRay) const
{
	path.clear();
	samplePath(path, startRay, 0);
}

vector<Path> MCRenderer::samplePathList(const vector<Ray>& startRayList) const
{
	Ray termRay;
	termRay.origin = vec3f(0, 0, 0);
	termRay.direction = vec3f(0, 0, 0);
	termRay.color = vec3f(0, 0, 0);
	termRay.directionSampleType = Ray::DEFINITE;
	termRay.directionProb = 1;
	termRay.insideObject = termRay.contactObject = termRay.intersectObject = NULL;

	vector<Path> pathList(startRayList.size());

	vector<Ray> rays;
	vector<unsigned> pathIDs;

	vector<Ray> candidateRays = startRayList;
	vector<unsigned> candidatePathIDs(startRayList.size());

	for(unsigned i=0; i<candidatePathIDs.size(); i++)
	{
		candidatePathIDs[i] = i;
	}

	for(int depth=0; depth<maxDepth && candidateRays.size()>0; depth++)
	{
		
		rays.clear();
		pathIDs.clear();
		for(int i=0; i<candidateRays.size(); i++)
		{
			pathList[candidatePathIDs[i]].push_back(candidateRays[i]);
			if(candidateRays[i].direction.length() > 0.5)
			{
				rays.push_back(candidateRays[i]);
				pathIDs.push_back(candidatePathIDs[i]);
			}
		}
		
		
		renderer->scene.fillIntersectObject(rays);

		
		candidateRays.clear();
		candidateRays.resize(rays.size());
#pragma omp parallel for
		for(int i=0; i<rays.size(); i++)
		{
			if(rays[i].insideObject && rays[i].intersectObject)
			{
				candidateRays[i] = rays[i].insideObject->scatter(rays[i]);
			}
			else if(rays[i].intersectObject)
			{
				candidateRays[i] = rays[i].intersectObject->scatter(rays[i]);
			}
			else
			{
				candidateRays[i] = termRay;
			}
			candidateRays[i].current_tid = rays[i].intersect_tid;
		}
		candidatePathIDs = pathIDs;

	}
	
	for (int i=0; i<rays.size(); i++)
	{
		if(rays[i].direction.length() > 0.5)
		{
			pathList[pathIDs[i]].push_back(termRay);
		}
	}
	
	return pathList;
}

void MCRenderer::showCurrentResult(const vector<vec3f>& pixelColors , unsigned* time)
{
	IplImage* image = cvCreateImage(cvSize(renderer->camera.width, renderer->camera.height), IPL_DEPTH_32F, 3);
	for(int x=0; x<renderer->camera.width; x++)
	{
		for(unsigned y=0; y<renderer->camera.height; y++)
		{
			vec3f rgb = pixelColors[y*renderer->camera.width + x];
			vec3f &bgr = ((vec3f*)image->imageData)[y*renderer->camera.width + x];
			bgr = vec3f(rgb.z, rgb.y, rgb.x);
		}
	}

	cvSetMouseCallback("Renderer", mouseEvent, this);
	if(savePath != "")
	{
		saveImagePFM(savePath, image);
	}

	if (time)
	{
		char timeStr[100];
		memset(timeStr , 0 , sizeof(timeStr));
		itoa(*time , timeStr , 10);
		string fileName("");
		for (int i = 0; i < savePath.length(); i++)
		{
			if (savePath[i] == '.')
				break;
			fileName.push_back(savePath[i]);
		}
		fileName.push_back('_');
		for (int i = 0; i < strlen(timeStr); i++)
			fileName.push_back(timeStr[i]);
		fileName += ".pfm";
		saveImagePFM(fileName , image);
	}

	for(int p=0; p<3*image->width*image->height; p++)
		((float*)image->imageData)[p] = powf(((float*)image->imageData)[p], 1/2.2);

	//cvShowImage("Renderer", image);
	
	//response(image);
	cvReleaseImage(&image);

	switch(showEvent)
	{
	case SHOW_PATH:
		if(showPath.size() > 1)
		{
			for(unsigned i=0; i<showPath.size()-1; i++)
			{
				vec3f point1 = showPath[i].origin;
				vec3f point2 = showPath[i+1].origin;
				vec2<float> pixel1 = renderer->getCamera().transToPixel(point1);
				vec2<float> pixel2 = renderer->getCamera().transToPixel(point2);
				CvPoint p1, p2;
				p1.x = pixel1.x;
				p1.y = pixel1.y;
				p2.x = pixel2.x;
				p2.y = pixel2.y;
				if(i==0)
					p1 = p2;
				cvLine(image, p1, p2, cvScalar(0, 1, 0));
				cvCircle(image, p2, 5, cvScalar(0, 0, 1), 2);
			}
		}
		break;
	/*case SHOW_PIXEL_VALUE:
		printf("%lf, %lf, %lf\n", pixelColors[pathPixelID].x, pixelColors[pathPixelID].y, pixelColors[pathPixelID].z);
		int x = pathPixelID % renderer->camera.width;
		int y = pathPixelID / renderer->camera.width;
		cvRectangle(image, cvPoint(x-1, y-1), cvPoint(x+1, y+1), cvScalar(0, 0, 1));
		break;*/
	}
}

string MCRenderer::response(const IplImage* currentImage)
{
	int key = cvWaitKey(10);
	
	switch(key)
	{
	case 's':
		if(currentImage)
		{
			IplImage *temp = cvCreateImage(cvSize(currentImage->width, currentImage->height), IPL_DEPTH_32F, 3);
			cvConvertScale(currentImage, temp, 1);
			saveImagePFM(savePath, temp);
			cvReleaseImage(&temp);
		}
		break;
	case 'q':
		return "quit";
	}
	return "";
}

void MCRenderer::saveImagePFM(const string& fileName, const IplImage* image)
{
	FILE* file;
	fopen_s(&file, fileName.c_str(), "wb");

	fprintf(file, "PF\n%d %d\n-1.000000\n", image->width, image->height);

	const float* data = (float*)image->imageData;
	for(int j=image->height-1; j>=0; j--)
	{
		for(unsigned i=0; i<image->width; i++)
		{
			fwrite(&data[3*(j*image->width+i)+2], sizeof(float), 1, file);
			fwrite(&data[3*(j*image->width+i)+1], sizeof(float), 1, file);
			fwrite(&data[3*(j*image->width+i)], sizeof(float), 1, file);
		}
	}
	fclose(file);
}

vector<vector<unsigned>> MCRenderer::testPathListVisibility(const vector<Path>& startPathList, const vector<Path>& endPathList)
{
	vector<vector<unsigned>> visibilityList(startPathList.size());
	vector<Ray> rays;

	unsigned lastPathIndex = 0;
	
	for(unsigned i=0; i<startPathList.size(); i++)
	{
		visibilityList[i].resize(ceil(startPathList[i].size()*endPathList[i].size()/32.0), 0);
		for(unsigned j=0; j<startPathList[i].size(); j++)
		{
			for(unsigned k=0; k<endPathList[i].size(); k++)
			{
				const Ray& startRay = startPathList[i][j];
				const Ray& endRay = endPathList[i][k];
				if(startRay.directionSampleType == Ray::RANDOM && endRay.directionSampleType == Ray::RANDOM)
				{
					Ray visRay = startRay;
					visRay.direction = endRay.origin - startRay.origin;
					visRay.intersectDist = visRay.direction.length();
					visRay.direction.normalize();
					rays.push_back(visRay);
				}
			}
		}
		if(rays.size() > 1024*1024 || i == startPathList.size()-1)
		{
			vector<bool> v = renderer->scene.testVisibility(rays);
			unsigned count = 0;
			for(unsigned pi=lastPathIndex; pi<=i; pi++)
			{
				unsigned subSize = startPathList[pi].size()*endPathList[pi].size();
				for(unsigned vi=0; vi<subSize; vi++)
				{
					const Ray& startRay = startPathList[pi][vi/endPathList[pi].size()];
					const Ray& endRay = endPathList[pi][vi%endPathList[pi].size()];
					if(startRay.directionSampleType == Ray::RANDOM && endRay.directionSampleType == Ray::RANDOM)
					{
						visibilityList[pi][vi/32] += unsigned(v[count] ? 1 : 0) << (vi%32);
						count++;
					}
				}
			}

			lastPathIndex = i+1;
			rays.clear();
		}
	}
	return visibilityList;
}

void MCRenderer::preprocessEmissionSampler(){ renderer->scene.preprocessEmissionSampler(); }

void MCRenderer::preprocessOtherSampler() { renderer->scene.preprocessOtherSampler(); }

Ray MCRenderer::genEmissiveSurfaceSample() const { return renderer->scene.genEmissionSample(); }

Ray MCRenderer::genOtherSurfaceSample() const { return renderer->scene.genOtherSample(); }