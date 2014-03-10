#include "StdAfx.h"
#include "Texture.h"
#include "ConfigManager.h"

void Texture::loadFile(const string& fileName)
{
	if(fileName.substr(fileName.size()-3, 3) == "pfm")
	{
		char c;
		FILE* file = fopen(fileName.c_str(), "rb");
		fscanf(file, "P");
		fscanf(file, "%c", &c);
		if(c == 'f')
			hasColor = false;
		else if(c == 'F')
			hasColor = true;

		fscanf(file, "%d %d", &width, &height);
		
		if(hasColor)
			colors.resize(width*height);
		else
			values.resize(width*height);

		float form;
		fscanf(file, "%f", &form);

		fread(&c, sizeof(char), 1, file);
		if(c == 0x0D)
			fread(&c, sizeof(char), 1, file);
		if(hasColor)
		{
			fread(colors.data(), sizeof(vec3f), colors.size(), file);
		}
		else
		{
			fread(values.data(), sizeof(float), values.size(), file);
		}
		fclose(file);
	}
}

vec3f Texture::getColor(const vec3f& coord) const
{
	vec3f transCoord = coordTransform * vec4f(coord, 1);
	float u = transCoord.x*width - 0.5;
	float v = transCoord.y*height - 0.5;

	int left = floor(u);
	int right = left + 1;
	int bottom = floor(v);
	int top = bottom + 1;

	u -= left;
	v -= bottom;

	left %= width;
	right %= width;
	top %= height;
	bottom %= height;

	if(left < 0) left += width;
	if(right < 0) right += width;
	if(top < 0) top += height;
	if(bottom < 0) bottom += height;

	if(hasColor)
	{
		vec3f lbColor = colors[bottom*width+left];
		vec3f rbColor = colors[bottom*width+right];
		vec3f ltColor = colors[top*width+left];
		vec3f rtColor = colors[top*width+right];
		vec3f color = vec3f(0, 0, 0);
		color += lbColor*(1-u)*(1-v);
		color += rbColor*u*(1-v);
		color += ltColor*(1-u)*v;
		color += rtColor*u*v;

		return vec3f(colorTransform * vec4f(color, 1));
	}
	else
	{
		float lbColor = values[bottom*width+left];
		float rbColor = values[bottom*width+right];
		float ltColor = values[top*width+left];
		float rtColor = values[top*width+left];
		float color = 0;
		color += lbColor*(1-u)*(1-v);
		color += rbColor*u*(1-v);
		color += ltColor*(1-u)*v;
		color += rtColor*u*v;
		return vec3f(colorTransform * vec4f(color, color, color, 1));
	}
}

vec3f Texture::getGrad(const vec3f& coord) const
{
	vec3f transCoord = coordTransform * vec4f(coord, 1);
	float u = transCoord.x*width - 0.5;
	float v = transCoord.y*height - 0.5;

	int left = floor(u);
	int right = left + 1;
	int bottom = floor(v);
	int top = bottom + 1;

	u -= left;
	v -= bottom;

	left %= width;
	right %= width;
	top %= height;
	bottom %= height;

	if(left < 0) left += width;
	if(right < 0) right += width;
	if(top < 0) top += height;
	if(bottom < 0) bottom += height;

	vec3f image_grad(0, 0, 0);

	if(hasColor)
	{
		float lbColor = colors[bottom*width+left].x;
		float rbColor = colors[bottom*width+right].x;
		float ltColor = colors[top*width+left].x;
		image_grad.x = rbColor - lbColor;
		image_grad.y = ltColor - lbColor;
	}
	else
	{
		float lbColor = values[bottom*width+left];
		float rbColor = values[bottom*width+right];
		float ltColor = values[top*width+left];
		image_grad.x = rbColor - lbColor;
		image_grad.y = ltColor - lbColor;
	}

	vec3f uv_grad;

	float var_n = image_grad.length();
	var_n = (colorTransform * vec4f(var_n, var_n, var_n, 0)).x;
	image_grad.normalize();
	uv_grad = vec3f(inverse(coordTransform)*vec4f(image_grad, 0));
	uv_grad.z = var_n;
	return uv_grad;
}

void Texture::loadTextureFromXML(const ConfigManager* cm, xml_node<>* nodeTex)
{
	loadFile(cm->getFullPath(nodeTex->first_node("path")->value()));
	if(nodeTex->first_node("coordTransform"))
		setCoordTransform(readMatrix(nodeTex->first_node("coordTransform")->value()));
	if(nodeTex->first_node("colorTransform"))
		setCoordTransform(readMatrix(nodeTex->first_node("colorTransform")->value()));
	
	if(nodeTex->first_node("colorScale"))
	{
		float scale = atof(nodeTex->first_node("colorScale")->value());
		setColorTransform(colorTransform*matrix4<float>(scale,0,0,0,0,scale,0,0,0,0,scale,0,0,0,0,1));
	}
}

void Texture::toGrayScale()
{
	hasColor = false;
	if(colors.size() == 0 || values.size() > 0)
	{
		colors.clear();
		return;
	}
	values.resize(colors.size());
	for(unsigned i=0; i<colors.size(); i++)
		values[i] = (colors[i].x + colors[i].y + colors[i].z) / 3;
	colors.clear();
}