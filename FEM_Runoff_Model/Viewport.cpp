#include "Viewport.hpp"

OffScreenBuffer viewportBuffer;

Vector2 delta;
float scale = 0.5f;
float margin = 0.05f;
float * nodesMeshVerts = NULL;
unsigned int * trianglesIndices = NULL;
int renderVertsCount = 0;
int renderTrisIndicesCount = 0;

std::unordered_map <int, Layer>  layers;
Shader triangleShader, pointShader, lineShader;
Shader testSuperTriShader;
float superTriMeshVerts[9];

void SetupMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount)
{
	glGenBuffers(1, &(targetGLData->vertexBufferObject));
	glGenVertexArrays(1, &(targetGLData->vertexArrayObject));
	glGenBuffers(1, &(targetGLData->vertexArrayElementObject));
	
	if (verticesCount < 3)
		return;

	UpdateMesh(targetGLData, mesh, verticesCount, indices, indexCount);
}

void UpdateMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount)
{
	glBindVertexArray(targetGLData->vertexArrayObject);
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexArrayObject);
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexBufferObject);
	if (mesh != NULL && verticesCount > 0)
		glBufferData(GL_ARRAY_BUFFER, verticesCount * sizeof(float), mesh, GL_STATIC_DRAW);
	
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, targetGLData->vertexArrayElementObject);

	if (indices != NULL && indexCount > 0)
	{
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexCount * sizeof(unsigned int), indices, GL_STATIC_DRAW);
	}
}

bool SetupShaders(Shader * targetGLData, char const * vertexShaderSource, char const * fragmentShaderSource)
{
	bool status = true;

	targetGLData->vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glShaderSource(targetGLData->vertexShader, 1, &vertexShaderSource, NULL);
	glCompileShader(targetGLData->vertexShader);

	int result;
	glGetShaderiv(targetGLData->vertexShader, GL_COMPILE_STATUS, &result);
	if (result == GL_FALSE)
	{
		std::cout << "failed to compile vertex shader" << std::endl;
		int length;
		glGetShaderiv(targetGLData->vertexShader, GL_INFO_LOG_LENGTH, &length);
		char * message = new char[length];
		glGetShaderInfoLog(targetGLData->vertexShader, length, &length, message);
		std::cout << message << std::endl;
		delete[] message;
		status = false;
	}

	targetGLData->fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(targetGLData->fragmentShader, 1, &fragmentShaderSource, NULL);
	glCompileShader(targetGLData->fragmentShader);

	glGetShaderiv(targetGLData->fragmentShader, GL_COMPILE_STATUS, &result);
	if (result == GL_FALSE)
	{
		std::cout << "failed to compile fragment shader" << std::endl;
		int length;
		glGetShaderiv(targetGLData->vertexShader, GL_INFO_LOG_LENGTH, &length);
		char * message = new char[length];
		glGetShaderInfoLog(targetGLData->vertexShader, length, &length, message);
		std::cout << message << std::endl;
		delete[] message;
		status = false;
	}

	targetGLData->program = glCreateProgram();
	glAttachShader(targetGLData->program, targetGLData->vertexShader);
	glAttachShader(targetGLData->program, targetGLData->fragmentShader);
	glLinkProgram(targetGLData->program);
	glValidateProgram(targetGLData->program);

	return status;
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
	//glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); //auto mipmaping //triggers enum error

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
	
	static int minDimension = Min(minViewportWidth, viewportHeight);
	glViewport(0, 0, minDimension, minDimension);

	glClearColor(viewportBGColour.x, viewportBGColour.y, viewportBGColour.z, viewportBGColour.w);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_DEPTH_TEST);

	if (nodes.size() > 2) //test supertriangle
	{
		glUseProgram(testSuperTriShader.program);
		glBindVertexArray(layers[10].meshData.vertexArrayObject);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, layers[10].meshData.vertexArrayElementObject);
		glDrawArrays(GL_TRIANGLES, 0, 3);
	}

	glBindVertexArray(layers[1].meshData.vertexArrayObject);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, layers[1].meshData.vertexArrayElementObject);
	
	//draw triangles
	glUseProgram(triangleShader.program);
	glDrawElements(GL_TRIANGLES, renderTrisIndicesCount, GL_UNSIGNED_INT, 0);
	
	//draw nodes
	glUseProgram(pointShader.program);
	glDrawArrays(GL_POINTS, 0, renderVertsCount);

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

	ImVec2 pos = ImGui::GetCursorScreenPos();
	ImGui::GetWindowDrawList()->AddImage((void*)viewportBuffer.texture, //texture to render (from an offscreen buffer)
		ImVec2(viewPortDimensions.positionX, viewPortDimensions.positionY), //NW Corner
		ImVec2(viewPortDimensions.positionX + viewPortDimensions.width, viewPortDimensions.positionY + viewPortDimensions.height), //SE corner
		ImVec2(0, 1), ImVec2(1, 0)); //UV
	ImGui::End();
}

bool InitViewport()
{
	std::cout << "Initializaing viewport\n";
	if (!SetupOffScreenBuffer(&viewportBuffer, minViewportWidth, viewportHeight))
		return FAILED_VIEWPORT_CREATE;

	layers.insert({ 1, Layer("MainMesh", MeshData(), 1) });
	//layers.insert({ 2, Layer("Triangles", MeshData(), 2) });
	//layers.insert({ 5, Layer("Test", MeshData(), 5) });
	
	GLErrorCheck();
	std::cout << "Initializaing layers\n";
	for (auto it = layers.begin(); it != layers.end(); ++it)
	{
		//renderVertsCount = 3;
		//SetupMesh(&viewportGLData, TestTriangle::vertices, sizeof(TestTriangle::vertices), NULL, 0); //Startup with test triangle
		SetupMesh(&(it->second.meshData), NULL, 0, NULL, 0); //startup empty
	}
	GLErrorCheck();
	std::cout << "Initializaing shaders\n";
	if ( !SetupShaders(&pointShader, CircleShader::vertex_shader_text, CircleShader::fragment_shader_text) ||
		!SetupShaders(&triangleShader, TestTriangle::vertex_shader_text, TestTriangle::fragment_shader_text) )
	{
		//return false;
		std::cout << "Error compiling shaders!" << std::endl;
	}
	
	GLErrorCheck();
	std::cout << "Finished initializaing viewpoty\n";


	TestSetupSuperTriangleRender();
	return true;
}

void UpdateViewport()
{
	UpdateOffScreenBuffer(&viewportBuffer, viewportWidth, viewportHeight);
	UpdateCoordinateConversionParameters();
	UpdateNodes();
	UpdateTriangles();
	
	UpdateMesh(&(layers[1].meshData), nodesMeshVerts, renderVertsCount, trianglesIndices, renderTrisIndicesCount);
	TestUpdateSuperTriangle();
}

void UpdateNodes()
{
	if (nodes.size() < 2)
		return;

	if (nodesMeshVerts != NULL)
		delete[] nodesMeshVerts; 

	renderVertsCount = nodes.size() * 3;
	nodesMeshVerts = new float[nodes.size() * 3];

	int counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		Vector2 relativePos = NormalizeCoordinates(*it);
		nodesMeshVerts[counter] = relativePos.x;
		nodesMeshVerts[counter + 1] = relativePos.y;
		nodesMeshVerts[counter + 2] = 0.0f;
		counter += 3;
	}
}

void UpdateTriangles()
{
	std::cout << "Attempting to update triangles with " << triangles.size() << " new tris\n";
	if (triangles.size() < 1)
		return;

	if (trianglesIndices != NULL)
		delete[] trianglesIndices;
	renderTrisIndicesCount = triangles.size() * 3;
	trianglesIndices = new unsigned int[renderTrisIndicesCount];

	int counter = 0;
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		trianglesIndices[counter] = it->second.vertIDs[0];
		trianglesIndices[counter + 1] = it->second.vertIDs[1];
		trianglesIndices[counter + 2] = it->second.vertIDs[2];
		counter += 3;
	}
}

void UpdateCoordinateConversionParameters()
{
	if (nodes.size() < 2)
		return;
	delta = nodesNE - nodesSW;
}

Vector2 NormalizeCoordinates(Vector2 & point)
{
	Vector2 normPos((2.0f - 2.0f * margin) * (point.x - nodesSW.x) / (delta.x) - 1.0f + margin,
					(2.0f - 2.0f * margin) * (point.y - nodesSW.y) / (delta.y) - 1.0f + margin);
	
	normPos = normPos * scale;
	return normPos;
}


static const char* vertex_shader_text =
		"#version 330\n"
		"layout(location = 0) in vec4 pos;\
		void main(void)\
		{\
			gl_Position = pos;\
		}";

static const char* fragment_shader_text =
"void main(void)\
		{\
			gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);\
		}";

void TestSetupSuperTriangleRender()
{
	SetupShaders(&testSuperTriShader, vertex_shader_text, fragment_shader_text);
	layers.insert({ 10, Layer("SuperTriangle", MeshData(), 10) });

	TestUpdateSuperTriangle();
}

void TestUpdateSuperTriangle()
{
	if (nodes.size() < 2)
		return;

	Vector2 diff = nodesNE - nodesSW;
	float boundHalfWidth = (diff.x / 2.0f) + SUPER_TRIANGLE_PADDING;
	float boundHeight = diff.y + 2.0f * SUPER_TRIANGLE_PADDING;
	Vector2 superVerts[3]{	Vector2(nodesSW.x - SUPER_TRIANGLE_PADDING + boundHalfWidth,
									nodesNE.y + SUPER_TRIANGLE_PADDING + boundHeight),
							Vector2(nodesSW.x - SUPER_TRIANGLE_PADDING - boundHalfWidth,
									nodesSW.y - SUPER_TRIANGLE_PADDING),
							Vector2(nodesNE.x + SUPER_TRIANGLE_PADDING + boundHalfWidth,
									nodesSW.y - SUPER_TRIANGLE_PADDING) };

	int counter = 0;
	for (int i = 0; i < 9; i += 3)
	{
		Vector2 relativePos = NormalizeCoordinates(superVerts[counter]);

		superTriMeshVerts[i] = relativePos.x;
		superTriMeshVerts[i + 1] = relativePos.y;
		superTriMeshVerts[i + 2] = 0.0f;

		counter++;
	}

	SetupMesh(&(layers[10].meshData), superTriMeshVerts, 9, NULL, 0);
	GLErrorCheck();
}