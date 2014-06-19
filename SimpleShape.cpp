#include "StdAfx.h"
#include "SimpleShape.h"

vec3f SimpleShape::getCenter() const
{
	return (minCoord+maxCoord)/2;
}

float SimpleShape::getDiagLen() const
{
	return (maxCoord - minCoord).length();
}

float SimpleShape::getTriangleArea(unsigned ti) const
{
	vec3f vps[3];
	for(unsigned i=0; i<3; i++)
	{
		vps[i] = getWorldVertexPosition(faceVertexIndexList[ti][i]);
	}
	vec3f l0 = vps[1] - vps[0];
	vec3f l1 = vps[2] - vps[1];
	float c = l0.cross(l1).length();
	return c/2;
}

vec3f SimpleShape::genRandTrianglePosition(unsigned ti) const
{
	vec3f vps[3];
	for(unsigned i=0; i<3; i++)
	{
		vps[i] = getWorldVertexPosition(faceVertexIndexList[ti][i]);
	}
	float u1 = RandGenerator::genFloat();
	float u2 = RandGenerator::genFloat();
	float su1 = sqrtf(u1);
	float u = 1.f - su1;
	float v = u2 * su1;
	vec3f res = vps[0] * u + vps[1] * v + vps[2] * (1.f - u - v);
	return res;
}

void SimpleShape::loadShape(const string &fileName, bool normalize, vector<SimpleShape*>* ss)
{
	bool split = (ss!=NULL);
	vector<SimpleShape*> &splitedShapes = *ss;
	int ret;
	vertexList.clear();
	faceVertexIndexList.clear();
	string suffix = fileName.substr(fileName.length()-4,4);
	if(suffix == ".shp")
	{
		FILE* file;
		fopen_s(&file, fileName.c_str(),"rb");
		int size;
		fread(&size,sizeof(int),1,file);
		vertexList.resize(size);
		fread(vertexList.data(),sizeof(float),size,file);
		fread(&size,sizeof(int),1,file);
		faceVertexIndexList.resize(size);
		fread(faceVertexIndexList.data(),sizeof(unsigned int),size,file);
		fclose(file);
	}
	else if(suffix == ".obj")
	{
		char line[BUFFERSIZE];
		char attrib[BUFFERSIZE];
		char parms[3][BUFFERSIZE];
		FILE* file;
		fopen_s(&file, fileName.c_str(),"r");

		unsigned current_fi = vertexList.size();
		unsigned current_fni = vertexNormalList.size();
		unsigned current_fti = vertexTexCoordList.size();

		SimpleShape* shape = NULL;
		
		while(fgets(line,BUFFERSIZE,file))
		{
			if(line[0] == '#')
				continue;
			int num = sscanf_s(line,"%s %s %s %s",attrib,BUFFERSIZE,parms[0],BUFFERSIZE,parms[1],BUFFERSIZE,parms[2],BUFFERSIZE);

			if(strcmp("g",attrib)==0 && split)
			{
				if(shape)
				{
					splitedShapes.push_back(shape);
				}
				shape = new SimpleShape;
				shape->name = parms[0];
				current_fi = vertexList.size();
				current_fni = vertexNormalList.size();
				current_fti = vertexTexCoordList.size();
			}

			if(num!=4)
				continue;
			if(strcmp("v",attrib)==0)
			{
				vec3f vert;
				sscanf_s(parms[0],"%f",&vert.x,sizeof(float));
				sscanf_s(parms[1],"%f",&vert.y,sizeof(float));
				sscanf_s(parms[2],"%f",&vert.z,sizeof(float));
				vertexList.push_back(vert);
				if(split && shape)
					shape->vertexList.push_back(vert);
			}
			if(strcmp("f",attrib)==0)
			{
				vec3ui tri, vnTri, tTri;
				bool has_n = false;
				bool has_t = false;

				ret = sscanf_s(parms[0],"%d/%d/%d",&tri.x,&tTri.x,&vnTri.x,sizeof(unsigned));
				ret = sscanf_s(parms[1],"%d/%d/%d",&tri.y,&tTri.y,&vnTri.y,sizeof(unsigned));
				ret = sscanf_s(parms[2],"%d/%d/%d",&tri.z,&tTri.z,&vnTri.z,sizeof(unsigned));

				if(ret==1)
				{
					ret = sscanf_s(parms[0],"%d//%d",&tri.x,&vnTri.x,sizeof(unsigned));
					ret = sscanf_s(parms[1],"%d//%d",&tri.y,&vnTri.y,sizeof(unsigned));
					ret = sscanf_s(parms[2],"%d//%d",&tri.z,&vnTri.z,sizeof(unsigned));
					has_n = ret == 2;

					if (ret == 1)
					{
						ret = sscanf_s(parms[0],"%d",&tri.x,sizeof(unsigned));
						ret = sscanf_s(parms[1],"%d",&tri.y,sizeof(unsigned));
						ret = sscanf_s(parms[2],"%d",&tri.z,sizeof(unsigned));
						has_n = false;
					}		
				}
				else
				{
					has_n = ret == 3;
					has_t = ret >= 2;
				}

				tTri -= vec3ui(1, 1, 1);
				tri -= vec3ui(1,1,1);
				vnTri -= vec3ui(1,1,1);
				faceVertexIndexList.push_back(tri);
				if(split && shape)
					shape->faceVertexIndexList.push_back(tri - vec3ui(current_fi, current_fi, current_fi));
				if(has_n)
				{
					faceVertexNormalIndexList.push_back(vnTri);
					if(split && shape)
						shape->faceVertexNormalIndexList.push_back(vnTri - vec3ui(current_fni, current_fni, current_fni));
				}
				if(has_t)
				{
					faceVertexTexCoordIndexList.push_back(tTri);
					if(split && shape)
						shape->faceVertexTexCoordIndexList.push_back(tTri - vec3ui(current_fti, current_fti, current_fti));
				}
			}
			if(strcmp("vn",attrib)==0)
			{
				vec3f vn;
				sscanf_s(parms[0],"%f",&vn.x,sizeof(float));
				sscanf_s(parms[1],"%f",&vn.y,sizeof(float));
				sscanf_s(parms[2],"%f",&vn.z,sizeof(float));
				vertexNormalList.push_back(vn);
				if(split && shape)
					shape->vertexNormalList.push_back(vn);
			}
			if(strcmp("vt",attrib)==0)
			{
				vec3f vt;
				sscanf_s(parms[0],"%f",&vt.x,sizeof(float));
				sscanf_s(parms[1],"%f",&vt.y,sizeof(float));
				sscanf_s(parms[2],"%f",&vt.z,sizeof(float));
				vertexTexCoordList.push_back(vt);
				if(split && shape)
					shape->vertexTexCoordList.push_back(vt);
			}
		}
		fclose(file);
		if(split)
		{
			if(shape)
			{
				splitedShapes.push_back(shape);
			}
		}
	}
	
	if(normalize)
	{
		matrix4<float> unitizeMat = this->unitize();
		for(unsigned i=0; split && i<splitedShapes.size(); i++)
		{
			splitedShapes[i]->setTransform(unitizeMat);
			splitedShapes[i]->applyTransform();
			splitedShapes[i]->setTransform(transform);
		}

		matrix4<float> changeHandness;
		changeHandness.set_scale(vec3f(-1.f , 1.f , 1.f));
		matrix4<float> trans = transform * unitizeMat;
		printf("=================\n");
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				printf("%.6f " , trans.element(i , j));
			printf("\n");
		}
		
	}
	getBoundingBox(minCoord, maxCoord);
}

void SimpleShape::applyTransform()
{
	matrix4<float> normalMat = transpose(inverse(transform));
	for(unsigned i=0; i<vertexList.size(); i++)
		vertexList[i] = vec3f(transform*vec4<float>(vertexList[i], 1));
	for(unsigned i=0; i<vertexNormalList.size(); i++)
		vertexNormalList[i] = vec3f(normalMat*vec4<float>(vertexNormalList[i], 0));
	transform.make_identity();
}

void SimpleShape::getBoundingBox(vec3f &minCoord, vec3f &maxCoord)
{
	if(vertexList.size()<1)
		return;
	minCoord.x=maxCoord.x=vertexList[0].x;
	minCoord.y=maxCoord.y=vertexList[0].y;
	minCoord.z=maxCoord.z=vertexList[0].z;
	for(unsigned i=1;i<getVertexNum();i++)
	{
		float &x = vertexList[i].x;
		float &y = vertexList[i].y;
		float &z = vertexList[i].z;
		minCoord.x = min(x,minCoord.x);
		maxCoord.x = max(x,maxCoord.x);
		minCoord.y = min(y,minCoord.y);
		maxCoord.y = max(y,maxCoord.y);
		minCoord.z = min(z,minCoord.z);
		maxCoord.z = max(z,maxCoord.z);
	}
}

matrix4<float> SimpleShape::unitize()
{
	vec3f minCoord, maxCoord;
	getBoundingBox(minCoord,maxCoord);
	vec3f center = (minCoord+maxCoord)/2;
	vec3f delta = maxCoord - minCoord;
	float maxLen = max(delta.x,delta.y);
	maxLen = max(maxLen,delta.z);
	for(unsigned i=0; i<getVertexNum(); i++)
	{
		vertexList[i] -= center;
		vertexList[i] /= maxLen/2;
	}
	//return matrix4<float>(2/maxLen, 0, 0, -2/maxLen*center.x, 0, 2/maxLen, 0, -2/maxLen*center.y, 0, 0, 2/maxLen, -2/maxLen*center.z, 0, 0, 0, 1);
	return matrix4<float>(2/maxLen, 0, 0, 0, 0, 2/maxLen, 0, 0, 0, 0, 2/maxLen, 0, -2/maxLen*center.x, -2/maxLen*center.y, -2/maxLen*center.z, 1);
}

vec3f SimpleShape::getTexCoord(unsigned fi, const vec3f& position) const
{
	if(!vertexTexCoordList.size())
		return vec3f(0, 0, 0);

	vec3f vps[3], vts[3];
	for(unsigned i=0; i<3; i++)
	{
		vps[i] = getWorldVertexPosition(faceVertexIndexList[fi][i]);
		vts[i] = vertexTexCoordList[faceVertexTexCoordIndexList[fi][i]];
	}

	vec3f b1 = vps[1] - vps[0];
	vec3f b2 = vps[2] - vps[0];

	vec3f v = position - vps[0];
	float d12 = b1.dot(b2);
	float l1 = b1.length();
	float l2 = b2.length();
	float u2 = (v.dot(b1)*d12 - v.dot(b2)*l1*l1) / (d12*d12 - l1*l1*l2*l2);
	float u1 = (v.dot(b1)-u2*d12)/(l1*l1);
	vec3f texCoord = (1-u1-u2)*vts[0] + u1*vts[1] + u2*vts[2];

	return texCoord;
}

vec3f SimpleShape::getWorldNormal(unsigned fi, const vec3f& position, bool flat) const
{
	vec3f vps[3], vns[3];
	for(unsigned i=0; i<3; i++)
	{
		if (fi >= faceVertexIndexList.size())
			printf("get world normal error , %d , %d\n" , fi , faceVertexIndexList.size());

		vps[i] = getWorldVertexPosition(faceVertexIndexList[fi][i]);
	}
	matrix4<float> normalMat = transpose(inverse(transform));
	vec3f b1 = vps[1] - vps[0];
	vec3f b2 = vps[2] - vps[0];
	vec3f faceNormal = b1.cross(b2);
	faceNormal.normalize();
	if(flat)
		return faceNormal;
	for(unsigned i=0; i<3; i++)
	{
		if(vertexNormalList.size())
		{
			vns[i] = vertexNormalList[faceVertexNormalIndexList[fi][i]];
			vns[i] = vec3f(normalMat * vec4<float>(vns[i], 0));
		}
		else
		{
			return faceNormal;
		}
	}
	
	vec3f v = position - vps[0];
	float d12 = b1.dot(b2);
	float l1 = b1.length();
	float l2 = b2.length();
	float u2 = (v.dot(b1)*d12 - v.dot(b2)*l1*l1) / (d12*d12 - l1*l1*l2*l2);
	float u1 = (v.dot(b1)-u2*d12)/(l1*l1);
	vec3f normal = (1-u1-u2)*vns[0] + u1*vns[1] + u2*vns[2];
	normal.normalize();
	return normal;
}

void SimpleShape::saveShape(const string &fileName)
{
	FILE* file;
	fopen_s(&file, fileName.c_str(),"wb");
	int size;

	size = vertexList.size()*3;
	fwrite(&size,sizeof(int),1,file);
	fwrite(vertexList.data(),sizeof(float),size,file);

	size = faceVertexIndexList.size()*3;
	fwrite(&size,sizeof(int),1,file);
	fwrite(faceVertexIndexList.data(),sizeof(unsigned int),size,file);

	fclose(file);
}

LocalFrame SimpleShape::getAutoGenWorldLocalFrame(unsigned fi, const vec3f& position, bool flat) const
{
	LocalFrame lf;
	lf.n = getWorldNormal(fi, position, flat);
	lf.buildFromNormal(lf.n);
	/*
	vec3f axis = vec3f(0, 1, 0).cross(lf.n);
	float angle = acos(clampf(vec3f(0, 1, 0).dot(lf.n), -1, 1));
	axis.normalize();
	lf.s = vec3f(rotMat(axis, angle)*vec4<float>(vec3f(1,0,0), 0));
	lf.t = vec3f(rotMat(axis, angle)*vec4<float>(vec3f(0,0,1), 0));
	*/
	return lf;
}