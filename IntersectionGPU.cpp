#include "stdafx.h"

#include "IntersectionGPU.h"

int IntersectionGPU::maxTexSize = 2048;
bool IntersectionGPU::renderToFBO = false;
int IntersectionGPU::width = 128;
int IntersectionGPU::height = 128;
GLuint IntersectionGPU::fboID = 0;
GLuint IntersectionGPU::color_rboID = 0;
GLuint IntersectionGPU::depth_rboID = 0;
Shader IntersectionGPU::shader;
vector<SimpleShape*> IntersectionGPU::shapes;
vector<IntersectionGPU::Ray> IntersectionGPU::rayList;

void IntersectionGPU::init()
{
	int ac = 0;
	char *av;
	glutInit(&ac, &av);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(width, height);
	glutCreateWindow("Renderer");
	glutDisplayFunc(render);
	/*glutMouseFunc(mouseFunc);
	glutMotionFunc(mouseActiveMotionFunc);
	glutMouseWheelFunc(mouseWheelFunc);
	glutKeyboardFunc(keyFunc);
	glutPassiveMotionFunc(mousePassiveMotionFunc);*/
	glutHideWindow();
	shader.init();
	shader.createProgram("Intersection", "Shader/Intersection.vert", "Shader/Intersection.frag");
	if(renderToFBO)
		initFBO();

	int size;
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxTexSize);
	glGetIntegerv(GL_MAX_VIEWPORT_DIMS, &size);
	maxTexSize = min2(maxTexSize, size);
	maxTexSize >>= 2;
	maxTexSize = min2(maxTexSize, 2048);

}

void IntersectionGPU::render()
{
	glViewport(0, 0, width, height);
	glutHideWindow();
	glClearColor(1, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	shader.useProgram("Intersection");

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, width, 0, height);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glBegin(GL_QUADS);
		glVertex2f(0, height);
		glVertex2f(width, height);
		glVertex2f(width, 0);
		glVertex2f(0, 0);
	glEnd();

	if(!renderToFBO)
		glutSwapBuffers();
	glFlush();
}

void IntersectionGPU::initFBO(int w, int h)
{
	if(w == 0)
		w = width;
	if(h == 0)
		h = height;

	renderToFBO = true;
	// Create
	// frame buffer object
	glGenFramebuffers(1, &fboID);
	glBindFramebuffer(GL_FRAMEBUFFER, fboID);

	// depth buffer with render buffer object
	/*glGenRenderbuffers(1, &depth_rboID);
	glBindRenderbuffer(GL_RENDERBUFFER, depth_rboID);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);*/

	// color buffer with texture object
	glGenTextures(1, &color_rboID);
	glBindTexture(GL_TEXTURE_2D, color_rboID);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, 0);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	//Attach depth buffer to FBO
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_rboID);

	// Attach color buffer to FBO
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_rboID, 0);
}

void IntersectionGPU::cleanUpFBO()
{
	renderToFBO = false;
	if (fboID)
	{
		glDeleteFramebuffers(1, &fboID);
		fboID = 0;
	}

	/*if (depth_rboID)
	{
		glDeleteRenderbuffers(1, &depth_rboID);
		depth_rboID = 0;
	})*/

	if (color_rboID)
	{
		glDeleteTextures(1, &color_rboID);
		color_rboID = 0;
	}
}

void IntersectionGPU::loadShapesToShader()
{
	vector<vec4f> vertexTex;
	vector<vec4f> triangleTex;
	vec4f offset(0, 0, 0, 0);
	for(unsigned i=0; i<shapes.size(); i++)
	{
		for(unsigned vi=0; vi<shapes[i]->vertexList.size(); vi++)
		{
			vec4f vertex;
			vec3f wvp = shapes[i]->getWorldVertexPosition(vi);
			vertex.x = wvp.x;
			vertex.y = wvp.y;
			vertex.z = wvp.z;
			vertexTex.push_back(vertex);
		}
		for(unsigned fi=0; fi<shapes[i]->faceVertexIndexList.size(); fi++)
		{
			vec4f tri;
			tri.x = shapes[i]->faceVertexIndexList[fi].x;
			tri.y = shapes[i]->faceVertexIndexList[fi].y;
			tri.z = shapes[i]->faceVertexIndexList[fi].z;
			triangleTex.push_back(tri+offset);
		}
		offset += shapes[i]->vertexList.size();
	}
	shader.useProgram("Intersection");
	shader.setUniform("nTriangles", int(triangleTex.size()));
	shader.setTexArray2D("vertexPositionTex", vertexTex, maxTexSize);
	shader.setTexArray2D("triangleVertexIndicesTex", triangleTex, maxTexSize);

	
}

vector<vec4f> IntersectionGPU::queryOneBatch(const vector<Ray>& rayList)
{
	vector<vec4f> originTex(rayList.size());
	vector<vec4f> directionTex(rayList.size());
	vector<vec4f> dist_tid(rayList.size());

	for(unsigned i=0; i<rayList.size(); i++)
	{
		originTex[i] = vec4f(rayList[i].origin, rayList[i].last_tid);
		directionTex[i] = vec4f(rayList[i].direction, 0);
	}

	shader.useProgram("Intersection");
	shader.setUniform("nQueries", int(rayList.size()));

	shader.setTexArray2D("queryRayOriginTex", originTex, maxTexSize);
	shader.setTexArray2D("queryRayDirectionTex", directionTex, maxTexSize);

	int prev_width = width;
	int prev_height = height;

	int w = rayList.size() < maxTexSize ? rayList.size() : maxTexSize;
	int h = ceil(rayList.size() / float(maxTexSize));

	width = w;
	height = h;

	initFBO();
	
	render();
	glReadBuffer(GL_COLOR_ATTACHMENT0);

	GLfloat* pixels = new GLfloat[4*w*h];
	glReadPixels(0, 0, w, h, GL_RGBA, GL_FLOAT, pixels);
	cleanUpFBO();
	
	shader.deleteLastTexture();
	shader.deleteLastTexture();

	memcpy(dist_tid.data(), pixels, rayList.size()*sizeof(vec4f));
	free(pixels);

	/*FILE* file = fopen("debug.txt", "w");
	for(unsigned i=0; i<dist_tid.size(); i++)
		fprintf(file, "%f, %f, %f, %f\n", dist_tid[i].x, dist_tid[i].y, dist_tid[i].z, dist_tid[i].w);
	fclose(file);*/

	width = prev_width;
	height = prev_height;

	return dist_tid;
}

vector<vec4f> IntersectionGPU::query(const vector<Ray>& rayList)
{

	vector<vec4f> dist_tid(rayList.size());
	int blockID = 0;
	int size = rayList.size();
	while(size > maxTexSize * maxTexSize)
	{
		vector<Ray> blockRays;
		vector<vec4f> blockResult;
		blockRays.assign(rayList.begin()+blockID*maxTexSize*maxTexSize, 
			rayList.begin()+(blockID+1)*maxTexSize*maxTexSize);
		blockResult = queryOneBatch(blockRays);
		copy(blockResult.begin(), blockResult.end(), dist_tid.begin()+blockID*maxTexSize*maxTexSize);
		blockID ++;
		size -= maxTexSize * maxTexSize;
	}
	
	vector<Ray> blockRays;
	vector<vec4f> blockResult;
	blockRays.assign(rayList.end() - size, rayList.end());
	blockResult = queryOneBatch(blockRays);
	copy(blockResult.begin(), blockResult.end(), dist_tid.end() - size);

	return dist_tid;
}

void IntersectionGPU::loadKDTreeToShader(const KDTree& tree)
{
	shader.useProgram("Intersection");

	vector<vec4f> nodes, nodes_minCoords, nodes_maxCoords, leaf_v1, leaf_v2, leaf_v3;
	tree.serializeForGPU(nodes, nodes_minCoords, nodes_maxCoords, leaf_v1, leaf_v2, leaf_v3);
	shader.setTexArray2D("kdTreeNodes_pi_ri_ls_le", nodes, maxTexSize);
	shader.setTexArray2D("kdTreeNodes_minCoords", nodes_minCoords, maxTexSize);
	shader.setTexArray2D("kdTreeNodes_maxCoords", nodes_maxCoords, maxTexSize);
	shader.setTexArray2D("leafVertexPosition1", leaf_v1, maxTexSize);
	shader.setTexArray2D("leafVertexPosition2", leaf_v2, maxTexSize);
	shader.setTexArray2D("leafVertexPosition3", leaf_v3, maxTexSize);
}