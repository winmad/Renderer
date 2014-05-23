#pragma once

#include "Renderer.h"
#include "macros.h"
#include <io.h>
#include <math.h>

typedef vector<Ray> Path;

class Renderer;

class MCRenderer
{
friend void mouseEvent (int evt, int x, int y, int flags, void* param);

public:

	enum {SHOW_PATH, SHOW_PIXEL_VALUE, NO_EVENT} showEvent;

	struct IntersectInfo
	{
		SceneObject *intersectObject;
		unsigned triangleID;
		float dist;
	};

	

protected:

	Renderer *renderer;

	void samplePath(Path& path, Ray& prevRay, unsigned depth, bool firstDiff = false) const;

	string savePath;

	vec4<float> connectColorProb(const Path& connectedPath, int connectIndex, bool merged = false);

protected:

	string response(const IplImage* currentImage = NULL);

	bool useMerge;

	unsigned maxDepth;

	unsigned pathPixelID;

	vector<Ray> showPath;

	void samplePath(Path& path, Ray& startRay) const;

	vector<Path> samplePathList(const vector<Ray>& startRayList) const;

	void showCurrentResult(const vector<vec3f>& pixelColors , unsigned* time = NULL);

	void eliminateVignetting(vector<vec3f>& pixelColors);

	SceneObject* findInsideObject(const vec3f& origin, const vec3f& direction);

	IntersectInfo intersect(const vec3f& origin, const vec3f& direction, KDTree::Condition* condition = NULL);

	void preprocessEmissionSampler();

	void preprocessOtherSampler();

	void preprocessVolumeSampler();

	Ray genEmissiveSurfaceSample(bool isUniform = false) const;

	Ray genOtherSurfaceSample(bool isUniform = false) const;

	Ray genVolumeSample(bool isUniform = false) const;

	void resetInsideObject(Ray& ray);

	bool testVisibility(const Ray& startRay, const Ray& endRay);

	vector<vector<unsigned>> testPathListVisibility(const vector<Path>& startPathList, const vector<Path>& endPathList);

	void saveImagePFM(const string& fileName, const IplImage* image);

	bool connectRays(Path& path, int connectIndex, bool merged = false);

public:

	void setMaxDepth(const int _maxDepth)
	{
		maxDepth = (unsigned)_maxDepth;
	}

	int getMaxDepth()
	{
		return maxDepth;
	}

	void setSavePath(const string& savePath);

	MCRenderer(Renderer *renderer)
	{
		this->renderer = renderer; 
		maxDepth = 20;
		pathPixelID = -1;
		useMerge = false;
	}

	virtual void preprocess(){}

	virtual vector<vec3f> renderPixels(const Camera& camera)
	{
		return vector<vec3f>(camera.generateRays().size(), vec3f(0, 0, 0));
	}

};

