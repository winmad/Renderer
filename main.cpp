// Renderer.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "SimpleShape.h"
#include "Renderer.h"

#include "IntersectionGPU.h"

int _tmain(int argc, _TCHAR* argv[])
{
	Renderer renderer;
	renderer.showWindow();
	renderer.waitForCommand();
	return 0;
}
