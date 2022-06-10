#pragma once
#include "GUI_Requirements.hpp"

extern OffScreenBuffer viewportBuffer;
extern std::list<OffScreenBuffer> offScreenBuffersList;

void SetupMesh(GLData * targetGLData, float * mesh, unsigned int verticesCount);
void UpdateMesh(GLData * targetGLData, float * mesh, unsigned int verticesCount);
void SetupShaders(GLData * targetGLData, char * vertexShaderSource, char * fragmentShaderSource);
bool SetupOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY); //resets active buffer to main buffer when done.
void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY); //resets active buffer to main buffer when done.

void RenderViewport(); //Renders viewport content to an offscreen buffer.
void DrawViewport(); //IMGUI painting commands for viewport window

void DeleteAllOffscreenBuffers();