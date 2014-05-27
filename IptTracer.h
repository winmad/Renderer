#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include <stack>
#include <utility>
#include "PointKDTree.h"
#include "CountHashGrid.h"
#include "macros.h"

struct IptPathState
{
	vec3f throughput , indirContrib;
	//vec3f contribs[4]; 
	Ray *ray , *lastRay , *originRay;
	int pathLen;
	vec3f pos;
	int index;
	//double accuProb;
};

class IptTracer : public MCRenderer
{
protected:
	unsigned spp;

	int mergeIterations;

	vector<vector<int> > partPathMergeIndex;

	vector<float> weights;

	int *revIndex;

	vector<IptPathState> partialSubPathList;

	vector<int> vis;

	stack<int> cycle;
	vector<pair<int , int> > edgeToRemove;

	bool dfs(int cur , int col);

	void movePaths(omp_lock_t& cmdLock , vector<Path>& , vector<Path*>&);

	void genLightPaths(omp_lock_t& cmdLock , vector<Path*>&);

	void genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>&);

	void mergePartialPaths(omp_lock_t& cmdLock);

	Ray genIntermediateSamples(Scene& scene);

	void mergePartialPaths(vector<vec3f>& contribs , const IptPathState& interState);

	vec3f colorByRayMarching(Path& eyeMergePath , PointKDTree<IptPathState>& partialSubPaths);

	vec3f colorByConnectingLights(Ray lastRay , Ray ray);

	void sampleMergePath(Path &path, Ray &prevRay, uint depth);

public:
	Real mergeRadius;
	Real lightMergeKernel , interMergeKernel;
	Real alpha;
	Real totArea , totVol;
	Real initialProb;
	unsigned timeInterval , lastTime;
	bool useWeight , usePPM , useDirIllu;

	int pixelNum , lightPathNum , cameraPathNum , interPathNum , partialPathNum;
	int lightPhotonNum , partialPhotonNum;

	IptTracer(Renderer* renderer) : MCRenderer(renderer)
	{ 
		alpha = 2.f / 3.f;
		spp = -1; 
		initialProb = 1.f;
		timeInterval = lastTime = 3600;

		pixelNum = renderer->camera.height * renderer->camera.width;
		cameraPathNum = pixelNum;
		lightPathNum = pixelNum / 2;
		interPathNum = pixelNum / 2;
		partialPathNum = interPathNum;

		usePPM = false;
		useDirIllu = false;
		if (usePPM)
		{
			mergeIterations = 0;
			useWeight = false;
		}
		else
		{
			mergeIterations = 1000;
			useWeight = true;
		}
	}
	void setRadius(const Real& r) { mergeRadius = r; }
	void setInitProb(const Real& r) { initialProb = r; }
	virtual vector<vec3f> renderPixels(const Camera& camera);
	
	Real connectFactor(Real pdf)
	{
		return pdf;
	}

	Real mergeFactor(Real *volScale = NULL , Real *initProb = NULL , Real *dirProb = NULL)
	{
		Real s = 1.0;
		if (volScale)
			s = *volScale;
		if (initProb)
			s *= *initProb;
		if (dirProb)
			s *= *dirProb;
		Real res = M_PI * mergeRadius * mergeRadius * partialPathNum * s;
		return res;
	}
	
	Real kernel(Real distSqr , Real radiusSqr)
	{
		Real s = max(0.f , 1.f - distSqr / radiusSqr);
		return 3.f * s * s / M_PI;
	}
};


struct GatherQuery
{
	vec3f color;
	IptTracer *tracer;
	IptPathState* cameraState;
	bool constKernel;
	int mergeNum;

	GatherQuery(IptTracer* tracer) { this->tracer = tracer; mergeNum = 0; constKernel = true; }

	void process(IptPathState& lightState)
	{
		Real volMergeScale = 1;
		if (cameraState->ray->insideObject && !cameraState->ray->contactObject)
			volMergeScale = 4.f / 3.f * tracer->mergeRadius;
		
		Real originProb = 1.f / tracer->totArea;
		if (cameraState->ray->insideObject && cameraState->ray->contactObject == NULL)
		{
			if (lightState.ray->insideObject == NULL || 
				lightState.ray->contactObject != NULL ||
				cameraState->ray->insideObject != lightState.ray->insideObject ||
				!cameraState->ray->insideObject->canMerge ||
				!lightState.ray->insideObject->canMerge)
			{
				return;
			}
			volMergeScale = 4.f / 3.f * tracer->mergeRadius;
			originProb = 1.f / tracer->totVol;
		}
		else if (cameraState->ray->contactObject)
		{
			if (lightState.ray->contactObject != cameraState->ray->contactObject ||
				!lightState.ray->contactObject->canMerge ||
				!cameraState->ray->contactObject->canMerge)
			{
				return;
			}
		}
		else 
		{
			return;
		}
		
		Ray outRay;
		vec3f bsdfFactor;
		
		outRay = *lightState.ray;
		outRay.direction = -(cameraState->lastRay->direction);
		bsdfFactor = lightState.lastRay->getBSDF(outRay);
		if (y(bsdfFactor) < 1e-7f)
			return;

		vec3f totContrib(0.f);
		if (lightState.index < tracer->lightPhotonNum)
			totContrib = lightState.throughput;
		else
			totContrib = lightState.indirContrib;
		//totContrib = lightState.indirContrib;

		vec3f tmp = totContrib * bsdfFactor * cameraState->throughput; 

		//Real lastPdf , weightFactor;
		//lastPdf = lightState.lastRay->getDirectionSampleProbDensity(outRay);
		//lastPdf = cameraState->accuProb;
		//originProb = cameraState->accuProb;
		//originProb = tracer->getOriginProb(tracer->countHashGrid , outRay.origin);

		//weightFactor = (tracer->mergeFactor(&volMergeScale , &originProb , &INV_2_PI)) /
		//	(tracer->connectFactor(lastPdf) + tracer->mergeFactor(&volMergeScale , &originProb , &INV_2_PI));

		//fprintf(fp , "weight = %.8f , cameraAccuProb = %.8f , hashOriginProb = %.8f\n" , weightFactor , cameraState->accuProb , originProb);

		//mergeNum++;

		vec3f res;
		float mergeKernel , pathNum;
		if (lightState.index < tracer->lightPhotonNum)
		{
			mergeKernel = tracer->lightMergeKernel / volMergeScale;
			pathNum = tracer->lightPathNum;
		}
		else
		{
			mergeKernel = tracer->interMergeKernel / volMergeScale;
			pathNum = tracer->partialPathNum;
		}

		if (constKernel)
		{	
			res = tmp * mergeKernel;
		}
		else
		{
			float distSqr = (cameraState->pos - lightState.pos).length();
			distSqr = distSqr * distSqr;
			mergeKernel = tracer->kernel(distSqr , tracer->mergeRadius * tracer->mergeRadius);
			if (abs(volMergeScale - 1.f) > 1e-6f)
				mergeKernel /= pathNum * tracer->mergeRadius * tracer->mergeRadius * tracer->mergeRadius;
			else 
				mergeKernel /= pathNum * tracer->mergeRadius * tracer->mergeRadius;

			res = tmp * mergeKernel;
		}
			
		//res *= weightFactor;

		if (!tracer->usePPM && tracer->useDirIllu && lightState.index < tracer->lightPhotonNum && 
			lightState.pathLen == 1 && lightState.ray->contactObject && !lightState.ray->contactObject->isVolumetric())
		{
			float dist = (lightState.lastRay->origin - lightState.ray->origin).length();
			float cosToLight = clampf(lightState.ray->getContactNormal().dot(-lightState.lastRay->direction) , 0.f , 1.f);
			float p1 = lightState.lastRay->originProb;
			float p2 = lightState.lastRay->originProb * lightState.lastRay->directionProb * cosToLight / (dist * dist) *
				M_PI * tracer->mergeRadius * tracer->mergeRadius * tracer->partialPathNum;
			float weightFactor = p2 / (p1 + p2);

			//printf("merge weight = %.8f\n" , weightFactor);

			res *= weightFactor;
		}

		color += res;
		/*
		if (cameraState->ray->insideObject)
		{
			fprintf(fp , "=====================\n");
			if (volMergeScale == 1)
				fprintf(fp , "surface\n");
			else 
				fprintf(fp , "volume\n");
			fprintf(fp , "res = (%.8f,%.8f,%.8f) \ntotContrib = (%.8f,%.8f,%.8f), bsdf = (%.8f,%.8f,%.8f),\n" ,
				res[0] , res[1] , res[2] ,
				totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2]);
		}
		*/
	}
};