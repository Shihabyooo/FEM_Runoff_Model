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
#include "FileIO.hpp"

#define MIN_VIEWPORT_DELTA 1.0

#define CLEAR_ARRAY(x) if (x != NULL) { delete[] x; } x = NULL;

enum ToolMode
{
	None = 0,
	ViewNode = 1,
	ViewElement = 2,
	MoveNode = 3
};

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

struct Icon
{
public:
	Icon()
	{

	}

	Icon(Image const & image)
	{
		glGenTextures(1, &texPtr);
		glBindTexture(GL_TEXTURE_2D, texPtr);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		if (image.bitmap == NULL)
			return;
		
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image.bitmap);
		width = image.width;
		height = image.height;
		isInit = true;
	}

	~Icon()
	{
		if (isInit)
			glDeleteTextures(0, &texPtr);
	}

	void operator= (Icon && icon2)
	{
		if (isInit)
			glDeleteTextures(0, &texPtr);


		texPtr = icon2.texPtr;
		width = icon2.width;
		height = icon2.height;

		isInit = icon2.isInit;
		icon2.isInit = false;
	}

	GLuint texPtr;
	int width = 0, height = 0;
	bool isInit = false;
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

extern ToolMode activeTool;

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