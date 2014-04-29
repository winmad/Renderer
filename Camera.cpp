#include "StdAfx.h"
#include "Camera.h"

Ray Camera::generateRay(unsigned pixelID) const
{
	vec3f front = focus - position;
	front.normalize();
	vec3f right = front.cross(up);
	right.normalize();
	vec3f orthoUp = right.cross(front);
	Ray ray;
	int x = pixelID % width;
	int y = pixelID / width;

	float u = RandGenerator::genFloat();
	float v = RandGenerator::genFloat();

	float xc = x + 0.5 - (width/2.0);
	float yc = (height/2.0) - 0.5 - y;

	ray.direction = front*sightDist + right * (x + u - (width/2.0)) + orthoUp * ((height/2.0) - v - y);
	ray.direction.normalize();
	ray.color = vec3f(1, 1, 1);
	ray.origin = position;
	ray.contactObject = (SceneObject*)this;

	if(!scene->usingGPU())
	{
		Scene::ObjSourceInformation osi;
		ray.intersectDist = scene->intersect(ray, osi);
		ray.intersectObjectTriangleID = osi.triangleID;
		ray.intersectObject = (ray.intersectDist >= 0) ? scene->objects[osi.objID] : NULL;
	}

	ray.insideObject = scene->findInsideObject(ray);
	float dist = sqrt(xc*xc+yc*yc+sightDist*sightDist);

	ray.directionProb = powf(dist, 3) / sightDist;

	ray.directionSampleType = Ray::RANDOM;
	ray.originSampleType = Ray::DEFINITE;
	ray.pixelID = pixelID;
	return ray;
}

float Camera::get_dw(unsigned pixelID) const
{
	int x = pixelID % width;
	int y = pixelID / width;
	float xc = x + 0.5 - (width/2.0);
	float yc = (height/2.0) - 0.5 - y;
	float dist = sqrt(xc*xc+yc*yc+sightDist*sightDist);
	return sightDist / powf(dist, 3);
}

vec3f Camera::eliminateVignetting(const vec3f& color, unsigned pixelID) const
{
	vec3f front = focus - position;
	front.normalize();
	vec3f right = front.cross(up);
	right.normalize();
	vec3f orthoUp = right.cross(front);
	Ray ray;
	int x = pixelID % width;
	int y = pixelID / width;

	float xc = x + 0.5 - (width/2.0);
	float yc = (height/2.0) - 0.5 - y;

	vec3f dir = front*sightDist + right * (x + 0.5 - (width/2.0)) + orthoUp * ((height/2.0) - 0.5 - y);
	dir.normalize();
	return color / powf(dir.dot(front), 4);
}


float Camera::getDirectionSampleProbDensity(const Ray& inRay, const Ray& outRay) const
{
	int x = outRay.pixelID % width;
	int y = outRay.pixelID / width;
	float xc = x + 0.5 - (width/2.0);
	float yc = (height/2.0) - 0.5 - y;
	float dist = sqrt(xc*xc+yc*yc+sightDist*sightDist);
	return powf(dist, 3) / sightDist;
}

vec3f Camera::getWorldNormal(unsigned fi, const vec3f& position, bool flat) const
{
	vec3f front = focus - position;
	front.normalize();
	return front;
}

vector<vec3f> Camera::generateRays() const
{
	vector<vec3f> rays(width*height);
	vec3f front = focus - position;
	front.normalize();
	vec3f right = front.cross(up);
	right.normalize();
	vec3f orthoUp = right.cross(front);
	for(unsigned y=0; y<height; y++)
	{
		for(unsigned x=0; x<width; x++)
		{
			rays[y*width+x] = front*sightDist + right * (x + 0.5 - (width/2.0)) + orthoUp * ((height/2.0) - 0.5 - y);
			rays[y*width+x].normalize();
		}
	}
	return rays;
}

vec2<float> Camera::transToPixel(const vec3f& point) const
{
	vec2<float> pixel;
	vec3f v = point - position;
	vec3f front = focus - position;
	front.normalize();
	vec3f right = front.cross(up);
	right.normalize();
	vec3f orthoUp = right.cross(front);
	float frontDist = v.dot(front);
	float ratio = sightDist / frontDist;
	pixel.x = v.dot(right) * ratio;
	pixel.y = v.dot(orthoUp) * ratio;
	pixel.x += width/2.0 - 0.5;
	pixel.y = height/2.0 - pixel.y - 0.5;
	return pixel;
}