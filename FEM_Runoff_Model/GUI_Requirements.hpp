#pragma once
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <glew.h>
#include <glfw3.h>
#include <list>
#include <iostream> //todo remove iostream and couts after implementing a decent logging functionality

//#include "Globals.hpp"
#include "ModelInterface.hpp"

struct Shader
{
public:
	GLuint program;
	GLuint vertexShader;
	GLuint fragmentShader;
};

struct MeshData //holds shaders, 
{
public:
	GLuint vertexBufferObject;
	GLuint vertexArrayObject;
	GLuint vertexArrayElementObject;
};

struct OffScreenBuffer
{
public:
	GLuint fbo;
	GLuint rbo;
	GLuint texture;
};

struct WindowDimensions
{
public:
	WindowDimensions()
	{

	};
	
	WindowDimensions(int _positionX, int _positionY, int _width, int _height)
	{
		positionX = _positionX;
		positionY = _positionY;
		width = _width;
		height = _height;
	};

	int positionX;
	int positionY;
	int width;
	int height;
};

void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY);
void RecomputeWindowElementsDimensions();// int newMainWinWidth, int newMainWinHeight);
void GLErrorCheck();

static const int minMainWinWidth = 800, minMainWinHeight = 600;
static const int minViewportWidth = 800, minViewportHeight = 600;

extern int mainWinWidth, mainWinHeight;
extern int viewportWidth, viewportHeight;

extern GLFWwindow * mainWindow;
//extern GLData viewportGLData;

extern ImVec4 mainBGColour;
extern ImVec4 viewportBGColour;

extern WindowDimensions leftPaneDimensions, logPaneDimensions, viewPortDimensions;

namespace CircleShader
{
	static const char* vertex_shader_text =
		"#version 330\n"
		"layout(location = 0) in vec4 pos;\
		void main(void)\
		{\
			gl_PointSize = 10.0;\
			gl_Position = pos;\
		}";

	static const char* fragment_shader_text =
		"#version 330\n"
		"void main(void)\
		{\
			gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);\
		}";
}

//For testing only:
namespace TestTriangle
{
	static float vertices[] = {
	0.0f,  0.5f,  0.0f,
	0.5f, -0.5f,  0.0f,
	-0.5f, -0.5f,  0.0f
	};

	static const char* vertex_shader_text =
		"#version 330\n"
		"in vec3 vp;"
		"void main() {"
		"  gl_Position = vec4(vp, 1.0);"
		"}";

	static const char* fragment_shader_text =
		"#version 330\n"
		"out vec4 frag_colour;"
		"void main() {"
		"  frag_colour = vec4(0.5, 0.0, 0.5, 1.0);"
		"}";
}