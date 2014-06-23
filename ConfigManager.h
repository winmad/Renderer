#pragma once
#include <io.h>
#include "textfile.h"
#include <rapidxml.hpp>
#include <string>
#include <unordered_map>
#include "xmlHelper.h"
#include "Renderer.h"

using namespace std;
using namespace rapidxml;

class Renderer;

class ConfigManager
{
private:
	string getFolder(const string& str) const
	{
		unsigned pos1 = str.rfind('/');
		unsigned pos2 = str.rfind('\\');
		unsigned pos = pos1 < pos2 ? pos1 : pos2;
		if(pos < str.size())
			return str.substr(0, pos);
		return "";
	}

	string rootFolder;

	string currentPath;

	unordered_map<string, pair<xml_document<>*, char*>> path_doc;

	Renderer* renderer;

	xml_node<>* findNode(const string& filePath, const string& nodeTag, const string& nodeName);

	xml_node<>* findNode(xml_node<>* root, const string& nodeTag, const string& nodeName);

	pair<string, string> getPathAndName(xml_node<>* node);

	void clear()
	{
		for(unordered_map<string, pair<xml_document<>*, char*>>::iterator it=path_doc.begin(); it!=path_doc.end(); it++)
		{
			delete (it->second).first;
			free(it->second.second);
		}
		path_doc.clear();
	}

	SceneObject* generateSceneObject(xml_node<>* nodeObj, xml_node<>* nodeMat, SimpleShape* shape = NULL);
public:
	string getFullPath(const string& fn) const
	{
		if(access(fn.c_str(), 0) == 0)
			return fn;
		string fullPath = getFolder(rootFolder + "/" + currentPath) + "/" + fn;
		if(access(fullPath.c_str(), 0) == 0)
			return fullPath;
		fullPath = rootFolder + "/" + fn;
		if(access(fullPath.c_str(), 0) == 0)
			return fullPath;
		return "";
	}
	ConfigManager(Renderer* renderer);
	void load(const string &configFilePath);
	string getRootPath() const { return rootFolder; }
};

