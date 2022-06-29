#pragma once
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <glew.h>
#include <glfw3.h>
#include <list>
#include <vector>
#include <string>>
#include <wchar.h>
#include <iostream> //todo remove iostream and couts after implementing a decent logging functionality

#include "ModelInterface.hpp"

#define MIN_VIEWPORT_DELTA 1.0

#define CLEAR_ARRAY(x) if (x != NULL) { delete[] x; } x = NULL;

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
		positionX = positionY = width = height = 0;
	};
	
	WindowDimensions(int _positionX, int _positionY, int _width, int _height)
	{
		positionX = _positionX;
		positionY = _positionY;
		width = _width;
		height = _height;
	};

	void SetDimensions(Vector2Int const & dimension)
	{
		width = dimension.x;
		height = dimension.y;
	}

	Vector2 Dimension()
	{
		return Vector2(width, height);
	}

	bool Contains(ImVec2 position)
	{
		return	position.x > positionX &&
				position.y > positionY &&
				position.x < positionX + width &&
				position.y < positionY + height;
	}

	Vector2 LocalPosFromGlobal(ImVec2 & const globalPos) //to convert positions sampled from dear imgui to ones relative to this subwindow
	{
		return Vector2(globalPos.x - positionX, globalPos.y - positionY);
	}

	int positionX;
	int positionY;
	int width;
	int height;
};

void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY);
void RecomputeWindowElementsDimensions();// int newMainWinWidth, int newMainWinHeight);
void GLErrorCheck();

static const int minMainWinWidth = 1024, minMainWinHeight = 768; //todo convert to Vector2Int
//static const int minViewportWidth = 800, minViewportHeight = 600; //todo convert to Vector2Int
//extern int mainWinWidth, mainWinHeight; //todo convert to Vector2Int

extern GLFWwindow * mainWindow;

extern ImVec4 mainBGColour;
extern ImVec4 viewportBGColour;

extern WindowDimensions leftPaneDimensions, logPaneDimensions, viewportDimensions;
extern WindowDimensions toolbarDimensions, statusBarDimensions;
extern Vector2D lastViewportSize;

namespace PointShader //TODO fix this
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
		"uniform vec4 diffuseCol;"
		"void main(void)\
		{\
			gl_FragColor = diffuseCol;\
			vec2 circularCoordinate = 2.0 * gl_PointCoord - 1.0;\
			if(dot(circularCoordinate, circularCoordinate) > 1.0)\
			{\
				discard;\
			}\
		}";
}

namespace PolygonShader
{
	static const char* vertex_shader_text =
		"#version 330\n"
		"in vec3 vp;"
		"void main() {"
		"  gl_Position = vec4(vp, 1.0);"
		"}";

	static const char* fragment_shader_text =
		"#version 330\n"
		"uniform vec4 diffuseCol;"
		"out vec4 frag_colour;"
		"void main() {"
		"  frag_colour = diffuseCol;"
		"}";
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