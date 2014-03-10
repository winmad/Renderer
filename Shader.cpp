#include "StdAfx.h"
#include "Shader.h"

void Shader::init()
{
	glewInit();
	if (glewIsSupported("GL_VERSION_2_0"))
		printf("Ready for OpenGL 2.0\n");
	else {
		printf("OpenGL 2.0 not supported\n");
		exit(1);
	}
	currentProgName = "";
	maxTexUnit = 0;
}

bool Shader::createProgram(const string& progName, const vector<string>& vertFileNameList, const vector<string>& fragFileNameList)
{
	bool allSuccess = true;
	if(progName == "")
		return false;
	vector<GLuint> vertShaderList(vertFileNameList.size());
	vector<GLuint> fragShaderList(fragFileNameList.size());
	GLuint progID = glCreateProgram();
	progName2progID[progName] = progID;
	for(unsigned i=0; i<vertShaderList.size(); i++)
	{
		vertShaderList[i] = glCreateShader(GL_VERTEX_SHADER);
		char *src = textFileRead(vertFileNameList[i].c_str());
		const char *src_const = src;
		glShaderSource(vertShaderList[i], 1, &src_const, NULL);
		free(src);
		glCompileShader(vertShaderList[i]);
		GLint success; 
		glGetShaderiv(vertShaderList[i], GL_COMPILE_STATUS, &success); 
		if (!success)
		{ 
			GLchar errorLog[1024]; 
			glGetShaderInfoLog(vertShaderList[i], sizeof(errorLog), NULL, errorLog); 
			printf("Error compiling %s: \n%s\n\n", vertFileNameList[i].c_str(), errorLog); 
			allSuccess = false;
			goto ret;
		}
		glAttachShader(progID, vertShaderList[i]);
		
	}
	for(unsigned i=0; i<fragShaderList.size(); i++)
	{
		fragShaderList[i] = glCreateShader(GL_FRAGMENT_SHADER);
		char *src = textFileRead(fragFileNameList[i].c_str());
		const char *src_const = src;
		glShaderSource(fragShaderList[i], 1, &src_const, NULL);
		free(src);
		glCompileShader(fragShaderList[i]);
		GLint success; 
		glGetShaderiv(fragShaderList[i], GL_COMPILE_STATUS, &success); 
		if (!success)
		{ 
			GLchar errorLog[1024]; 
			glGetShaderInfoLog(fragShaderList[i], sizeof(errorLog), NULL, errorLog); 
			printf("Error compiling %s:\n%s\n\n", fragFileNameList[i].c_str(), errorLog); 
			allSuccess = false;
			goto ret;
		}
		glAttachShader(progID, fragShaderList[i]);
		
	}
	glLinkProgram(progID);
	GLint linkSuccess;
	
	glGetProgramiv(progID, GL_LINK_STATUS, &linkSuccess); 
	
	if (!linkSuccess)
	{ 
		GLchar errorLog[1024]; 
		glGetProgramInfoLog(progID, sizeof(errorLog), NULL, errorLog); 
		printf("Error linking shader program %s:\n%s\n\n", progName.c_str(), errorLog); 
		
		goto ret;
	}
	
ret:
	
	for(unsigned i=0; i<vertShaderList.size(); i++)
	{
		if(vertShaderList[i])
			glDeleteShader(vertShaderList[i]);
	}
	for(unsigned i=0; i<fragShaderList.size(); i++)
	{
		if(fragShaderList[i])
			glDeleteShader(fragShaderList[i]);
	}
	if(!allSuccess)
		deleteProgram(progName);
	else
		currentProgName = progName;
	return allSuccess;
}

bool Shader::useProgram(const string& progName)
{
	if(progName == "")
		glUseProgram(0);
	else
		glUseProgram(progName2progID[progName]);
	currentProgName = progName;
	return true;
}

bool Shader::deleteProgram(const string& progName)
{
	GLuint &progID = progName2progID[progName];
	if(progID)
		glDeleteProgram(progID);
	progName2progID.erase(progName);
	return true;
}

void Shader::setUniform(const string& progName, const string& varName, const int v)
{
	glUniform1i(glGetUniformLocation(progName2progID[progName], varName.c_str()), v);
}

void Shader::setUniform(const string& progName, const string& varName, const float v)
{
	glUniform1f(glGetUniformLocation(progName2progID[progName], varName.c_str()), v);
}

void Shader::setUniform(const string& progName, const string& varName, const vec2f& v)
{
	glUniform2f(glGetUniformLocation(progName2progID[progName], varName.c_str()), v.x, v.y);
}

void Shader::setUniform(const string& progName, const string& varName, const vec3f& v)
{
	glUniform3f(glGetUniformLocation(progName2progID[progName], varName.c_str()), v.x, v.y, v.z);
}

void Shader::setUniform(const string& progName, const string& varName, const vec4f& v)
{
	glUniform4f(glGetUniformLocation(progName2progID[progName], varName.c_str()), v.x, v.y, v.z, v.w);
}

void Shader::setUniform(const string& progName, const string& varName, const matrix4<float>& v)
{
	glUniformMatrix4fv(glGetUniformLocation(progName2progID[progName], varName.c_str()), 1, false, (GLfloat*)&v);
}

Shader::~Shader(void)
{
	for(unordered_map<string, GLuint>::iterator it=progName2progID.begin(); it!=progName2progID.end(); it++)
		glDeleteProgram(it->second);
}

void Shader::setTexArray1D(const string& varName, const vector<vec4f>& tex)
{
	GLuint texID = 0;
	glActiveTexture(GL_TEXTURE0+maxTexUnit);
	glGenTextures(1, &texID);
	texIDs.push_back(texID);
	glBindTexture(GL_TEXTURE_1D, texID);

	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA32F, tex.size(), 0, GL_RGBA, GL_FLOAT, tex.data());

	setUniform(varName, int(maxTexUnit));
	maxTexUnit++;

	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glActiveTexture(GL_TEXTURE0+maxTexUnit);
}

void Shader::setTexArray2D(const string& varName, const vector<vec4f>& tex, int width)
{
	if(tex.size() < width)
		width = tex.size();

	GLuint texID = 0;
	glActiveTexture(GL_TEXTURE0+maxTexUnit);
	glGenTextures(1, &texID);
	texIDs.push_back(texID);
	glBindTexture(GL_TEXTURE_2D, texID);

	int height = ceil(tex.size() / float(width));

	GLfloat *pixels = new GLfloat[width*height*4];
	memcpy(pixels, tex.data(), tex.size()*sizeof(vec4f));
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, pixels);
	free(pixels);

	setUniform(varName, int(maxTexUnit));
	maxTexUnit++;

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glActiveTexture(GL_TEXTURE0+maxTexUnit);
}

void Shader::setTexArray2D(const string& varName, const vector<vector<vec4f>>& tex)
{
	GLfloat *pixels;
	GLuint texID = 0;
	glActiveTexture(GL_TEXTURE0+maxTexUnit);
	glGenTextures(1, &texID);
	texIDs.push_back(texID);
	glBindTexture(GL_TEXTURE_2D, texID);

	unsigned height = tex.size();
	unsigned width = tex[0].size();

	pixels = new GLfloat[width*height*4];
	for(unsigned i=0; i<height; i++)
	{
		memcpy(pixels+i*width*4, tex[i].data(), width*sizeof(vec4f));
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, pixels);
	free(pixels);

	setUniform(varName, int(maxTexUnit));
	maxTexUnit++;

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glActiveTexture(GL_TEXTURE0+maxTexUnit);
}