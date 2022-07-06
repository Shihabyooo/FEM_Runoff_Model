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
extern Vector2D viewBounds[2];

void DrawViewport(); //IMGUI painting commands for viewport window
bool InitViewport();

void UpdateViewport();
void SetViewBounds(Vector2D swCorner, Vector2D nwCorner);
