#include "Viewport.hpp"

OffScreenBuffer viewportBuffer;
std::list<OffScreenBuffer> offScreenBuffersList;

void SetupMesh(GLData * targetGLData, float * mesh, unsigned int verticesCount)
{
	glGenBuffers(1, &(targetGLData->vertexBufferObject));
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexBufferObject);
	glBufferData(GL_ARRAY_BUFFER, verticesCount * sizeof(float), mesh, GL_STATIC_DRAW);

	glGenVertexArrays(1, &(targetGLData->vertexArrayObject));
	glBindVertexArray(targetGLData->vertexArrayObject);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexArrayObject);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void UpdateMesh(GLData * targetGLData, float * mesh, unsigned int verticesCount)
{
	std::cout << "Attemptint to update Mesh\n";
	//glGenBuffers(1, &(targetGLData->vertexBufferObject));
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexBufferObject);
	glBufferData(GL_ARRAY_BUFFER, verticesCount * sizeof(float), mesh, GL_STATIC_DRAW);

	//glGenVertexArrays(1, &(targetGLData->vertexArrayObject));
	glBindVertexArray(targetGLData->vertexArrayObject);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexArrayObject);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void SetupShaders(GLData * targetGLData, char * vertexShaderSource, char * fragmentShaderSource)
{
	targetGLData->vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(targetGLData->vertexShader, 1, &vertexShaderSource, NULL);
	glCompileShader(targetGLData->vertexShader);

	targetGLData->fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(targetGLData->fragmentShader, 1, &fragmentShaderSource, NULL);
	glCompileShader(targetGLData->fragmentShader);

	targetGLData->program = glCreateProgram();
	glAttachShader(targetGLData->program, targetGLData->vertexShader);
	glAttachShader(targetGLData->program, targetGLData->fragmentShader);
	glLinkProgram(targetGLData->program);
}

bool SetupOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY) //resets active buffer to main buffer when done.
{
	//custom framebuffer (texture type) to use for sub-window
	glGenFramebuffers(1, &(buffer->fbo));
	glBindFramebuffer(GL_FRAMEBUFFER, buffer->fbo); //bind framebuffer so following commands apply to it
	glDrawBuffer(GL_COLOR_ATTACHMENT0);

	glGenTextures(1, &(buffer->texture));
	glBindTexture(GL_TEXTURE_2D, buffer->texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, sizeX, sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); //auto mipmaping


	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, buffer->texture, 0); //bind to current active fbo
	glBindTexture(GL_TEXTURE_2D, 0); //return to default texture.

	glGenRenderbuffers(1, &(buffer->rbo));
	glBindRenderbuffer(GL_RENDERBUFFER, buffer->rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, sizeX, sizeY);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, buffer->fbo);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, buffer->rbo);

	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
	{
		std::cout << "Error creating custom framebuffer" << std::endl;
		return false;
	}

	offScreenBuffersList.push_back(*buffer);

	glBindFramebuffer(GL_FRAMEBUFFER, 0); //restore default framebuffer

	return true;
}

void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY) //resets active buffer to main buffer when done.
{
	glBindFramebuffer(GL_FRAMEBUFFER, buffer->fbo);
	glBindTexture(GL_TEXTURE_2D, buffer->texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, sizeX, sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glBindRenderbuffer(GL_RENDERBUFFER, buffer->rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, sizeX, sizeY);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void RenderViewport() //Renders viewport content to an offscreen buffer.
{
	//Draw to offscreen buffer.
	/*nodesMeshVerts = new float[18]{ 0.0f,1.0f,0.0f,	1.0f, 0.0f, 0.0f,	0.0f, 0.0f,0.0f,
									0.0f, 1.0f, 0.0f,	1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f };
	UpdateMesh(&viewportGLData, nodesMeshVerts, 18);
	renderVertsCount = 16;*/

	glBindFramebuffer(GL_FRAMEBUFFER, viewportBuffer.fbo);
	int minDimension = viewportWidth > viewportHeight ? viewportHeight : viewportWidth;
	glViewport(0, 0, minDimension, minDimension);
	glClearColor(viewportBGColour.x, viewportBGColour.y, viewportBGColour.z, viewportBGColour.w);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_DEPTH_TEST);

	glUseProgram(viewportGLData.program);
	glBindVertexArray(viewportGLData.vertexArrayObject);
	glDrawArrays(GL_TRIANGLES, 0, renderVertsCount);

	//restore screen buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void DrawViewport() //IMGUI painting commands for viewport window
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	ImGui::SetNextWindowPos(ImVec2(viewPortDimensions.positionX, viewPortDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(viewPortDimensions.width, viewPortDimensions.height), ImGuiCond_Always);
	//ImGui::SetNextWindowSizeConstraints(ImVec2(viewportWidth, viewportHeight), ImVec2(INFINITY, INFINITY));

	ImGui::Begin("Viewport", NULL, windowFlags);

	//ImVec2 currentViewportSize = ImGui::GetWindowSize();
	//viewportWidth = currentViewportSize.x < minViewportWidth ? minViewportWidth : currentViewportSize.x;
	//viewportHeight = currentViewportSize.y < minViewportHeight ? minViewportHeight : currentViewportSize.y;

	ImVec2 pos = ImGui::GetCursorScreenPos();
	ImGui::GetWindowDrawList()->AddImage((void*)viewportBuffer.texture, //texture to render (from an offscreen buffer)
		ImVec2(viewPortDimensions.positionX, viewPortDimensions.positionY), //NW Corner
		ImVec2(viewPortDimensions.positionX + viewPortDimensions.width, viewPortDimensions.positionY + viewPortDimensions.height), //SE corner
		ImVec2(0, 1), ImVec2(1, 0)); //UV
	ImGui::End();
}

//void DeleteOffscreenBuffer(GLuint * fbo);
void DeleteAllOffscreenBuffers()
{
	while (offScreenBuffersList.size() > 0)
	{
		glDeleteBuffers(1, &(offScreenBuffersList.back().fbo));
		glDeleteBuffers(1, &(offScreenBuffersList.back().rbo));
		offScreenBuffersList.pop_back();
	}
}