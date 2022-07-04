#pragma once
#include "GUI_Requirements.hpp"

//TODO replace unnecessarily exposed declarations here with forward declarations in source file.

struct Layer
{
public:
	Layer()
	{

	};

	Layer(std::string & _name, MeshData _data, int _order)
	{
		name = _name;
		meshData = _data;
		order = _order;
	};

	Layer(const char * _name, MeshData _data, int _order)
	{
		name = std::string(_name);
		meshData = _data;
		order = _order;
	};

	std::string name = "Layer";
	MeshData meshData;
	int order;
};

extern OffScreenBuffer viewportBuffer;

extern std::unordered_map <int, Layer> layers;
extern Shader triangleShader, pointShader, lineShader;

extern Vector2D currenViewportHoverPos; //in worldspace coordinates
extern Vector2D currenViewportHoverPosPixels; //in local pixel screenspace of viewport
extern bool isHoveringViewport;
extern double scale;
extern double screenAspectRatio;
extern double worldAspectRatio;
extern Vector2D viewBounds[2];

void SetupMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount); //test test
void UpdateMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount);
bool SetupShaders(Shader * targetGLData, char const * vertexShaderSource, char const * fragmentShaderSource);
bool SetupOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY); //resets active buffer to main buffer when done.
void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY); //resets active buffer to main buffer when done.

//void RenderViewport(); //Renders viewport content to an offscreen buffer.
void DrawViewport(); //IMGUI painting commands for viewport window

bool InitViewport();

void UpdateViewport();
void UpdateContent();

void UpdateNodes();
void UpdateTriangles();

void UpdateViewBounds();
void SetViewBounds(Vector2D swCorner, Vector2D nwCorner);
void UpdateCoordinateProjectionParameters();
Vector2D NormalizeCoordinates(Vector2D & point); //converts from world space to viewport space
Vector2D ScreenToWorldSpace(Vector2D screenCoordinates);
void PanView(Vector2D posDelta);


void TestSetupSuperTriangleRender();
void TestUpdateSuperTriangle();
