#include "StdAfx.h"
#include "ConfigManager.h"
#include "PathTracer.h"
#include "SceneDiffuseObject.h"
#include "SceneEmissiveObject.h"
#include "SceneRefractiveObject.h"
#include "SceneVPMObject.h"
#include "SceneHeterogeneousVolume.h"
#include "SceneReflectiveObject.h"
#include "BidirectionalPathTracer.h"
#include "SceneGlossyObject.h"
#include "SceneParallelLight.h"
#include "ScenePhongObject.h"
#include "Photonmap.h"
#include "VCMTracer.h"
#include "IptTracer.h"
#include "PathTracerTest.h"


ConfigManager::ConfigManager(Renderer* renderer)
{
	this->renderer = renderer;
}

SceneObject* ConfigManager::generateSceneObject(xml_node<>* nodeObj, xml_node<>* nodeMat, SimpleShape* shape)
{
	SceneObject* sceneObject = NULL;
	if(strcmp(nodeMat->first_node("type")->value(), "Diffuse") == 0)
	{
		SceneDiffuseObject *obj = new SceneDiffuseObject(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;
		obj->getMaterial()->loadMaterialFromXML(this, nodeMat);
	}
	if(strcmp(nodeMat->first_node("type")->value(), "Glossy") == 0)
	{
		SceneGlossyObject *obj = new SceneGlossyObject(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		obj->getMaterial()->loadMaterialFromXML(this, nodeMat);
		if(shape)
			*((SimpleShape*)obj) = *shape;
	}
	if(strcmp(nodeMat->first_node("type")->value(), "Phong") == 0)
	{
		ScenePhongObject *obj = new ScenePhongObject(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		obj->getGlossyMaterial()->loadMaterialFromXML(this, nodeMat->first_node("GlossyComponent"));
		obj->getDiffuseMaterial()->loadMaterialFromXML(this, nodeMat->first_node("DiffuseComponent"));
		if(shape)
			*((SimpleShape*)obj) = *shape;
	}
	if(strcmp(nodeMat->first_node("type")->value(), "VPM") == 0)
	{
		SceneVPMObject *obj = new SceneVPMObject(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		obj->setSigma(readVec(nodeMat->first_node("sigma_s")->value()),readVec(nodeMat->first_node("sigma_a")->value()));
		//obj->setSSS(std::string(nodeMat->first_node("sss")->value()));
		obj->setInRate(atof(nodeMat->first_node("inrate")->value()));
		obj->setG(atof(nodeMat->first_node("g")->value()));
		if (nodeMat->first_node("stepsize"))
			obj->setStepSize(atof(nodeMat->first_node("stepsize")->value()));
		if (nodeMat->first_node("scale"))
			obj->setScale(atof(nodeMat->first_node("scale")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;
	}

	if(strcmp(nodeMat->first_node("type")->value(), "GridVolume") == 0)
	{
		HeterogeneousVolume *obj = new HeterogeneousVolume(&renderer->scene, atof(nodeMat->first_node("g")->value()));
		sceneObject = obj;
		if(nodeMat->first_node("IOR"))
			obj->setIOR(atof(nodeMat->first_node("IOR")->value()));
		if(nodeMat->first_node("stepsize"))
			obj->setStepSize(atof(nodeMat->first_node("stepsize")->value()));
		if(nodeMat->first_node("sigma_s") && nodeMat->first_node("sigma_a")){
			obj->setSigma(readVec(nodeMat->first_node("sigma_s")->value()),readVec(nodeMat->first_node("sigma_a")->value()));
		}
		if(nodeMat->first_node("density_scale"))
			obj->setDensityScale(atof(nodeMat->first_node("density_scale")->value()));
		if(nodeMat->first_node("scattering_scale"))
			obj->setScatteringScale(atof(nodeMat->first_node("scattering_scale")->value()));
		if(nodeMat->first_node("absorption_scale"))
			obj->setAbsorptionScale(atof(nodeMat->first_node("absorption_scale")->value()));

		if(nodeMat->first_node("lookup_flag"))
			obj->setVolumeLookUp(atof(nodeMat->first_node("lookup_flag")->value()));


		if(nodeMat->first_node("subsurface"))
			obj->setSubSurface(atof(nodeMat->first_node("subsurface")->value()));

		if(nodeObj->first_node("transform"))
			obj->transform = (readMatrix(nodeObj->first_node("transform")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;

		std::cout << "GridVolume Obj " << std::endl;
		return obj;
		
	}

	if(strcmp(nodeMat->first_node("type")->value(), "Emissive") == 0)
	{
		SceneEmissiveObject *obj = new SceneEmissiveObject(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		obj->setColor(readVec(nodeMat->first_node("color")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;
	}
	if(strcmp(nodeMat->first_node("type")->value(), "ParallelLight") == 0)
	{
		SceneParallelLight *obj = new SceneParallelLight(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		obj->setColor(readVec(nodeMat->first_node("color")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;
	}
	if(strcmp(nodeMat->first_node("type")->value(), "Refractive") == 0)
	{
		SceneRefractiveObject *obj = new SceneRefractiveObject(&renderer->scene);
		sceneObject = obj;
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;
		if(nodeMat->first_node("coeff"))
			obj->setRefrCoeff(atof(nodeMat->first_node("coeff")->value()));
		if(nodeMat->first_node("surfColor"))
			obj->setSurfColor(readVec(nodeMat->first_node("surfColor")->value()));
		if(nodeMat->first_node("decayColor"))
			obj->setDecayColor(readVec(nodeMat->first_node("decayColor")->value()));
	}
	if(strcmp(nodeMat->first_node("type")->value(), "Reflective") == 0)
	{
		SceneReflectiveObject *obj = new SceneReflectiveObject(&renderer->scene);
		sceneObject = obj;
		if(nodeMat->first_node("color"))
			obj->setColor(readVec(nodeMat->first_node("color")->value()));
		if(nodeObj->first_node("transform"))
			obj->setTransform(readMatrix(nodeObj->first_node("transform")->value()));
		if(shape)
			*((SimpleShape*)obj) = *shape;
	}
	if(nodeMat->first_node("BumpTexture"))
	{
		sceneObject->bumpTex.loadTextureFromXML(this, nodeMat->first_node("BumpTexture"));
		sceneObject->bumpTex.toGrayScale();
	}
	return sceneObject;
}

void ConfigManager::load(const string &configFilePath)
{
	clear();
	
	pair<string, string> path_name;

	rootFolder = getFolder(configFilePath);

	xml_node<> *nodeConfig = findNode(configFilePath, "Config", "");
	
	if (nodeConfig->first_node("useGPU"))
	{
		int flag = atoi(nodeConfig->first_node("useGPU")->value());
		if (flag == 1)
			renderer->scene.setGPU(true);
	}

	xml_node<> *nodeSceneConfig = nodeConfig->first_node("Scene");

	path_name = getPathAndName(nodeSceneConfig);
	xml_node<> *nodeScene = findNode(path_name.first, "Scene", path_name.second);
	if(!nodeScene)
		nodeScene = findNode(configFilePath, "Scene", path_name.second);

	currentPath = path_name.first;

	for(xml_node<>* nodeGroupObj = nodeScene->first_node("GroupObject"); nodeGroupObj; nodeGroupObj = nodeGroupObj->next_sibling("GroupObject"))
	{
		string path = nodeGroupObj->first_node("filePath")->value();
		
		if(_access(path.c_str(), 0))
		{
			path = rootFolder + '/' + path;
		}

		SimpleShape ss;
		if(nodeGroupObj->first_node("transform"))
			ss.setTransform(readMatrix(nodeGroupObj->first_node("transform")->value()));
		vector<SimpleShape*> shapes;
		ss.loadShape(path, true, &shapes);

		unordered_map<string, SimpleShape*> name_shape;

		for(unsigned i=0; i<shapes.size(); i++)
			name_shape[shapes[i]->name] = shapes[i];
		
		for(xml_node<>* nodeObj = nodeGroupObj->first_node("Object"); nodeObj; nodeObj = nodeObj->next_sibling("Object"))
		{
			path_name = getPathAndName(nodeObj->first_node("Material"));
			xml_node<>* nodeMat = findNode(path_name.first, "Material", path_name.second);
			if(!nodeMat && nodeObj->first_node("Material")->first_node("name"))
				nodeMat = findNode(nodeScene, "Material", path_name.second);
			if(!nodeMat)
				nodeMat = nodeObj->first_node("Material");
			SimpleShape *shape = name_shape[nodeObj->first_attribute("Name")->value()];

			renderer->scene.objects.push_back(generateSceneObject(nodeObj, nodeMat, shape));
		}
		for(unsigned i=0; i<shapes.size(); i++)
			delete shapes[i];
	}

	for(xml_node<>* nodeObj = nodeScene->first_node("Object"); nodeObj; nodeObj = nodeObj->next_sibling("Object"))
	{
		path_name = getPathAndName(nodeObj->first_node("Material"));
		xml_node<>* nodeMat = findNode(path_name.first, "Material", path_name.second);
		if(!nodeMat && nodeObj->first_node("Material")->first_node("name"))
			nodeMat = findNode(nodeScene, "Material", path_name.second);
		if(!nodeMat)
			nodeMat = nodeObj->first_node("Material");

		string path = nodeObj->first_node("filePath")->value();
		if(_access(path.c_str(), 0))
		{
			path = rootFolder + '/' + path;
		}

		SceneObject *obj = generateSceneObject(nodeObj, nodeMat);

		obj->loadShape(path , true); // NEED TO BE TRUE WHEN NOT TESTING
		// FALSE WHEN test_cornell_box

		if(obj->isVolumetric()){
			if(nodeMat->first_node("density_file"))
				((HeterogeneousVolume*)obj)->loadDensityMap(nodeMat->first_node("density_file")->value());

			if(nodeMat->first_node("volume_data_s") && nodeMat->first_node("volume_data_a")){
				((HeterogeneousVolume*)obj)->loadSubSurfaceVolumeData(nodeMat->first_node("volume_data_s")->value(), nodeMat->first_node("volume_data_a")->value());
			}
			std::cout << "BBox = " << (vec3f(obj->transform * vec4f(obj->minCoord, 1))) << " " << (vec3f(obj->transform * vec4f(obj->maxCoord, 1))) << std::endl;
		}


		renderer->scene.objects.push_back(obj);
	}

	xml_node<> *nodeCam = findNode(nodeScene, "Camera", nodeSceneConfig->first_node("camera")->value());
	renderer->camera.width = atoi(nodeCam->first_node("width")->value());
	renderer->camera.height = atoi(nodeCam->first_node("height")->value());
	renderer->camera.sightDist = atoi(nodeCam->first_node("sightDist")->value());
	renderer->camera.focus = readVec(nodeCam->first_node("focus")->value());
	renderer->camera.position = readVec(nodeCam->first_node("position")->value());
	renderer->camera.up = readVec(nodeCam->first_node("up")->value());

	renderer->scene.buildKDTree();

	if(renderer->mcRenderer)
	{
		delete renderer->mcRenderer;
		renderer->mcRenderer = NULL;
	}

	if(nodeConfig->first_node("renderer"))
	{
		string typeName = nodeConfig->first_node("renderer")->value();
		if(typeName == "PathTracer" || typeName == "PT")
		{
			renderer->mcRenderer = new PathTracer(renderer);
		}
		if(typeName == "PTtest")
		{
			renderer->mcRenderer = new PathTracerTest(renderer);
		}
		if(typeName == "PhotonMap" || typeName == "PM")
		{
			renderer->mcRenderer = new PhotonMap(renderer);
			//static_cast<PhotonMap*>(renderer->mcRenderer)->setPhotonsWant(atoi(nodeConfig->first_node("photons")->value()));
			static_cast<PhotonMap*>(renderer->mcRenderer)->setRadius(atof(nodeConfig->first_node("radius")->value()));
		}
		if(typeName == "BidirectionalPathTracer" || typeName == "BPT")
		{
			renderer->mcRenderer = new BidirectionalPathTracer(renderer);
		}
		/*if (typeName == "newBPT")
		{
			renderer->mcRenderer = new NewBidirectionalPathTracer(renderer);
		}*/
		if(typeName == "VCMTracer" || typeName == "VCM")
		{
			renderer->mcRenderer = new VCMTracer(renderer);
			((VCMTracer*)renderer->mcRenderer)->setRadius(atof(nodeConfig->first_node("radius")->value()));
		}
		if(typeName == "IptTracer" || typeName == "IPT")
		{
			renderer->mcRenderer = new IptTracer(renderer);
			((IptTracer*)renderer->mcRenderer)->setRadius(atof(nodeConfig->first_node("radius")->value()));
			if (nodeConfig->first_node("gatherRadius"))
				((IptTracer*)renderer->mcRenderer)->gatherRadius = (atof(nodeConfig->first_node("gatherRadius")->value()));
		}
		if (nodeConfig->first_node("maxDepth"))
			renderer->mcRenderer->setMaxDepth(atoi(nodeConfig->first_node("maxDepth")->value()));
	}

	xml_node<>* nodeRenderer = nodeConfig->first_node("Renderer");

	if(nodeRenderer)
	{
		string typeName = nodeRenderer->first_node("type")->value();
		if(typeName == "PathTracer" || typeName == "PT")
		{
			renderer->mcRenderer = new PathTracer(renderer);
		}
		if(typeName == "PTtest")
		{
			renderer->mcRenderer = new PathTracerTest(renderer);
		}
		if(typeName == "PhotonMap" || typeName == "PM")
		{
			renderer->mcRenderer = new PhotonMap(renderer);
			//if(nodeRenderer->first_node("photons"))
			//	static_cast<PhotonMap*>(renderer->mcRenderer)->setPhotonsWant(atoi(nodeRenderer->first_node("photons")->value()));
			if(nodeRenderer->first_node("radius"))
				static_cast<PhotonMap*>(renderer->mcRenderer)->setRadius(atof(nodeRenderer->first_node("radius")->value()));
		}
		if(typeName == "BidirectionalPathTracer" || typeName == "BPT")
		{
			renderer->mcRenderer = new BidirectionalPathTracer(renderer);
		}
		//if (typeName == "newBPT")
		//{
		//	renderer->mcRenderer = new NewBidirectionalPathTracer(renderer);
		//}
		if(typeName == "VCMTracer" || typeName == "VCM")
		{
			renderer->mcRenderer = new VCMTracer(renderer);
			if(nodeRenderer->first_node("radius"))
				((VCMTracer*)renderer->mcRenderer)->setRadius(atof(nodeRenderer->first_node("radius")->value()));
		}
		if(typeName == "IptTracer" || typeName == "IPT")
		{
			renderer->mcRenderer = new IptTracer(renderer);
			if(nodeRenderer->first_node("radius"))
				((IptTracer*)renderer->mcRenderer)->setRadius(atof(nodeRenderer->first_node("radius")->value()));
			if (nodeRenderer->first_node("gatherRadius"))
				((IptTracer*)renderer->mcRenderer)->gatherRadius = (atof(nodeRenderer->first_node("gatherRadius")->value()));
			if (nodeRenderer->first_node("pathRatio"))
				((IptTracer*)renderer->mcRenderer)->pathRatio = (atof(nodeRenderer->first_node("pathRatio")->value())); 
			if (nodeRenderer->first_node("mergeRatio"))
				((IptTracer*)renderer->mcRenderer)->mergeRatio = (atof(nodeRenderer->first_node("mergeRatio")->value())); 
			if (nodeRenderer->first_node("usePPM"))
				((IptTracer*)renderer->mcRenderer)->usePPM = (atoi(nodeRenderer->first_node("usePPM")->value()) == 1 ? true : false);
			if (nodeRenderer->first_node("pathNum"))
				((IptTracer*)renderer->mcRenderer)->totPathNum = (atoi(nodeRenderer->first_node("pathNum")->value()));
			if (nodeRenderer->first_node("isDebug"))
				((IptTracer*)renderer->mcRenderer)->isDebug = (atoi(nodeRenderer->first_node("isDebug")->value()) == 1 ? true : false);
			if (nodeRenderer->first_node("useUniformSur"))
				((IptTracer*)renderer->mcRenderer)->useUniformSur = (atoi(nodeRenderer->first_node("useUniformSur")->value()) == 1 ? true : false);
			if (nodeRenderer->first_node("useDirIllu"))
				((IptTracer*)renderer->mcRenderer)->useDirIllu = (atoi(nodeRenderer->first_node("useDirIllu")->value()) == 1 ? true : false);
		}
		if (nodeRenderer->first_node("maxDepth"))
			renderer->mcRenderer->setMaxDepth(atoi(nodeRenderer->first_node("maxDepth")->value()));
	}

	if(nodeConfig->first_node("savePath"))
	{
		renderer->mcRenderer->setSavePath(nodeConfig->first_node("savePath")->value());
	}
}

xml_node<>* ConfigManager::findNode(const string& filePath, const string& nodeTag, const string& nodeName)
{
	if(filePath == "")
		return NULL;
	string fullFilePath = filePath;
	if(_access(filePath.c_str(), 0))
		fullFilePath = rootFolder + string("/") + filePath;
	xml_document<> *doc;
	if(path_doc.find(fullFilePath) != path_doc.end())
	{
		doc = path_doc[fullFilePath].first;
	}
	else
	{
		doc = new xml_document<>;
		char* text = textFileRead(fullFilePath.c_str());
		path_doc[fullFilePath] = make_pair(doc, text);
		doc->parse<0>(text);
	}
	for(xml_node<> *node=doc->first_node(nodeTag.c_str()); node; node = node->next_sibling(nodeTag.c_str()))
	{
		if(node->first_attribute("Name") == NULL || nodeName == "")
			return node;
		if(node->first_attribute("Name") && node->first_attribute("Name")->value() == nodeName)
			return node;
	}
	return NULL;
}

xml_node<>* ConfigManager::findNode(xml_node<>* root, const string& nodeTag, const string& nodeName)
{
	for(xml_node<> *node=root->first_node(nodeTag.c_str()); node; node = node->next_sibling(nodeTag.c_str()))
	{
		if(node->first_attribute("Name") == NULL || nodeName == "")
			return node;
		if(node->first_attribute("Name") && node->first_attribute("Name")->value() == nodeName)
			return node;
	}
	return NULL;
}

pair<string, string> ConfigManager::getPathAndName(xml_node<>* node)
{
	pair<string, string> path_name;
	path_name.first = node->first_node("path") ? node->first_node("path")->value() : currentPath;
	path_name.second = node->first_node("name") ? node->first_node("name")->value() : "";
	return path_name;
}