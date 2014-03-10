#pragma once
#include "nvVector.h"
#include "nvMatrix.h"
#include "xmlHelper.h"
#include <string>
#include <vector>
#include "EXT/rapidxml/rapidxml.hpp"

class ConfigManager;
using namespace rapidxml;
using namespace std;
using namespace nv;

class Texture
{
private:
	bool hasColor;
	int width, height;
	vector<float> values;
	vector<vec3f> colors;
	matrix4<float> coordTransform, colorTransform;
public:
	void setCoordTransform(const matrix4<float>& transform){ coordTransform = transform; }
	void setColorTransform(const matrix4<float>& transform){ colorTransform = transform; }
	void loadFile(const string& fileName);
	void setColor(const vec3f& color)
	{ 
		hasColor = true;
		colors.resize(1, color); 
		width = height = 1;
	}

	void toGrayScale();

	int size() const { return width * height; }

	vec3f getColor(const vec3f& coord) const;

	vec3f getGrad(const vec3f& coord) const;

	void loadTextureFromXML(const ConfigManager* cm, xml_node<>* nodeTex);
};

