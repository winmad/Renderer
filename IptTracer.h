#pragma once
#include "mcrenderer.h"
#include <omp.h>
#include "PointKDTree.h"
#include "CountHashGrid.h"
#include "macros.h"

struct IptPathState
{
	vec3f throughput , indirContrib;
	//vec3f contribs[4]; 
	Ray *ray , *lastRay , *originRay;
	bool isSpecularPath;
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

	CountHashGrid countHashGrid;

	void movePaths(omp_lock_t& cmdLock , vector<Path>& , vector<Path*>&);

	void genLightPaths(omp_lock_t& cmdLock , vector<Path*>&);

	void genIntermediatePaths(omp_lock_t& cmdLock , vector<Path*>&);

	void mergePartialPaths(omp_lock_t& cmdLock);

	Ray genIntermediateSamples(vector<IptPathState>& partialSubPathList , Scene& scene , int *index = NULL);

	Ray link(const Path& path, int i, int j);

	void calcEyePathProbs(Path& eyePath , vector<double>& probDir , vector<double>& probRev);

	void mergePartialPaths(vector<vec3f>& contribs , const IptPathState& lightState , const int mergeIters);

	vec3f colorByMergingPaths(IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths);

	vec3f colorByMergingSurface(IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths);

	vec3f colorByMergingVolume(IptPathState& cameraState, PointKDTree<IptPathState>& partialSubPaths , vec3f& tr);

	vec3f colorByRayMarching(Path& eyeMergePath , PointKDTree<IptPathState>& partialSubPaths);

	vec3f colorByConnectingLights(const Camera& camera, const IptPathState& cameraState);

	vec3f colorByConnectingCamera(const Camera& camera, const IptPathState& lightState , int& _x , int& _y);

	void sampleMergePath(Path &path, Ray &prevRay, uint depth);

public:
	Real mergeRadius;
	Real mergeKernel;
	Real alpha;
	Real totArea , totVol;
	Real initialProb;
	unsigned timeInterval , lastTime;
	bool useWeight , usePPM;

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
		lightPathNum = pixelNum;
		interPathNum = pixelNum;
		partialPathNum = pixelNum;

		usePPM = true;
		if (usePPM)
		{
			mergeIterations = 0;
			useWeight = false;
		}
		else
		{
			mergeIterations = 5;
			useWeight = true;
		}
	}
	void setRadius(const Real& r) { mergeRadius = r; }
	void setInitProb(const Real& r) { initialProb = r; }
	virtual vector<vec3f> renderPixels(const Camera& camera);
	
	Real getOriginProb(CountHashGrid& hashGrid , vec3f& pos , const bool isVol)
	{
		CountQuery query(pos);
		hashGrid.count(query);
		return query.count * (hashGrid.mInvCellSize * hashGrid.mInvCellSize * hashGrid.mInvCellSize);
	}
	
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
			volMergeScale = 4.0 / 3.0 * tracer->mergeRadius;
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
		totContrib = lightState.throughput + lightState.indirContrib;
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
		if (constKernel)
		{
			res = tmp * (tracer->mergeKernel / volMergeScale);
		}
		else
		{
			float distSqr = (cameraState->pos - lightState.pos).length();
			distSqr = distSqr * distSqr;
			float mergeKernel = tracer->kernel(distSqr , tracer->mergeRadius * tracer->mergeRadius);
			if (abs(volMergeScale - 1.f) > 1e-6f)
				mergeKernel /= (float)tracer->partialPathNum * tracer->mergeRadius * tracer->mergeRadius * tracer->mergeRadius;
			else 
				mergeKernel /= (float)tracer->partialPathNum * tracer->mergeRadius * tracer->mergeRadius;

			res = tmp * mergeKernel;
		}
			
		//res *= weightFactor;

		color += res;
		/*
		vec3f resx = tracer->renderer->camera.eliminateVignetting(res , cameraState->index) *
			tracer->pixelNum;	
		if (y(res) > 1)
		{
			fprintf(fp , "=====================\n");
			if (volMergeScale == 1)
				fprintf(fp , "surface\n");
			else 
				fprintf(fp , "volume\n");
			fprintf(fp , "res = (%.8f,%.8f,%.8f) \ntotContrib = (%.8f,%.8f,%.8f), bsdf = (%.8f,%.8f,%.8f),\n cameraThr = (%.8f,%.8f,%.8f) \nweightFactor = %.8f, originProb = %.8f, lastpdf = %.8f\n" ,
				res[0] , res[1] , res[2] ,
				totContrib[0] , totContrib[1] , totContrib[2] , bsdfFactor[0] , bsdfFactor[1] , bsdfFactor[2] , 
				cameraState->throughput[0] , cameraState->throughput[1] , cameraState->throughput[2] , 
				weightFactor , originProb , lastPdf);
		}
		*/
	}
};