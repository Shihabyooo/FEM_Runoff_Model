#include "Viewport.hpp"

OffScreenBuffer viewportBuffer;

Vector2 delta;
double scale = 0.5f;
double aspectRatio = 1.0f;
//float margin = 0.05f;
float * nodesMeshVerts = NULL;
unsigned int * trianglesIndices = NULL;
int renderVertsCount = 0;
int renderTrisIndicesCount = 0;

std::unordered_map <int, Layer>  layers;
Shader triangleShader, pointShader, lineShader;
Shader testSuperTriShader;
float superTriMeshVerts[9];

Vector2 lastViewportSize;
Vector2 viewBounds[2];

//Create vertex buffer (VBO), vertex array object (VAO) and vertex array element object, store id in targetGLData struct. \
if verticesCount is greater than or equal to 3, UpdateMesh() is called.
void SetupMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount)
{
	glGenBuffers(1, &(targetGLData->vertexBufferObject));
	glGenVertexArrays(1, &(targetGLData->vertexArrayObject));
	glGenBuffers(1, &(targetGLData->vertexArrayElementObject));
	
	if (verticesCount >= 3)
		UpdateMesh(targetGLData, mesh, verticesCount, indices, indexCount);
}

//Updates mesh objects with IDs set in targetGLData with mesh[] and indices[]. verticesCount must be <= #of floats in mesh. \
indexCount must <= # of ints in indices. This function does nothing if mesh = NULL or verticesCount < 1. \
if indices = NULL || indexCount < 1, no data will be set in the element array buffer.
void UpdateMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount)
{
	if (mesh == NULL || verticesCount < 1)
		return;
	glBindVertexArray(targetGLData->vertexArrayObject);
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexArrayObject);
	glBindBuffer(GL_ARRAY_BUFFER, targetGLData->vertexBufferObject);
	glBufferData(GL_ARRAY_BUFFER, verticesCount * sizeof(float), mesh, GL_STATIC_DRAW);
	
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, targetGLData->vertexArrayElementObject);

	if (indices != NULL && indexCount > 0)
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexCount * sizeof(unsigned int), indices, GL_STATIC_DRAW);
}

//Compiles vertex and fragment shaders, attaches them to an openGL program, and stores the IDs of all three in targetGLData. \
Returns false if compilation of either shader fails and posts message on the log if this occured.
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
		LogMan::Log("failed to compile vertex shader", LOG_ERROR);

		int length;
		glGetShaderiv(targetGLData->vertexShader, GL_INFO_LOG_LENGTH, &length);
		char * message = new char[length];
		glGetShaderInfoLog(targetGLData->vertexShader, length, &length, message);
		LogMan::Log(message, LOG_ERROR);
		delete[] message;
		status = false;
	}

	targetGLData->fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(targetGLData->fragmentShader, 1, &fragmentShaderSource, NULL);
	glCompileShader(targetGLData->fragmentShader);

	glGetShaderiv(targetGLData->fragmentShader, GL_COMPILE_STATUS, &result);
	if (result == GL_FALSE)
	{
		LogMan::Log("failed to compile fragment shader", LOG_ERROR);
		int length;
		glGetShaderiv(targetGLData->vertexShader, GL_INFO_LOG_LENGTH, &length);
		char * message = new char[length];
		glGetShaderInfoLog(targetGLData->vertexShader, length, &length, message);
		LogMan::Log(message, LOG_ERROR);
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

//Initializes an off-screen frame buffer object (fbo), creates a 2D RGB texture and a render buffer \
object (rbo) and attaches both to the fbo. Returns false if initialization fails. Resets bound framebuffer \
to default (0) before returning.
bool SetupOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY) 
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
		LogMan::Log("Failed to creat offscreen framebuffer", LOG_ERROR);
		return false;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0); //restore default framebuffer

	return true;
}

//Changes the dimensions of the texture assigned to the fbo to sizeX * sizeY. Resets active buffer to main buffer when done.
void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY)
{
	//Bind the target fbo, bind the texture of the fbo, create a new data for the texture with target res, do the same for rbo.
	glBindFramebuffer(GL_FRAMEBUFFER, buffer->fbo);
	glBindTexture(GL_TEXTURE_2D, buffer->texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, sizeX, sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glBindRenderbuffer(GL_RENDERBUFFER, buffer->rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, sizeX, sizeY);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

float scaleChangeTisk = 0.05f; //TODO make this value dynamically set based on current bounds.
void HandleMousePan()
{
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 * posAtClick = io.MouseClickedPos;
	if (ImGui::IsMouseDragging(2) && viewPortDimensions.Contains(posAtClick[2]))
	{
		ImVec2 mouseDelta = io.MouseDelta;
		PanView(Vector2(-1.0f * mouseDelta.x, mouseDelta.y) * scale);
	}
}

void HandleMouseZoom()
{
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 posAtScroll = io.MousePos;
	float scroll = io.MouseWheel;
	if (scroll != 0.0f && viewPortDimensions.Contains(posAtScroll))
	{
		std::cout << "scroll: " << scroll << std::endl;

		//posAtScroll is the target point
		//recompute viewBounds that maintain this point at its place
		//recompute projection parameters.
	}
}


//Renders the viewport content to the offscreen buffer.
void RenderViewport() //Renders viewport content to an offscreen buffer.
{
	HandleMousePan();
	HandleMouseZoom();
	////Test mesh (quad covering whole viewport).
	/*nodesMeshVerts = new float[18]{ 0.0f,1.0f,0.0f,	1.0f, 0.0f, 0.0f,	0.0f, 0.0f,0.0f,
									0.0f, 1.0f, 0.0f,	1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f };
	UpdateMesh(&viewportGLData, nodesMeshVerts, 18);
	renderVertsCount = 16;*/

	//Objects in the back are drawn first.

	glBindFramebuffer(GL_FRAMEBUFFER, viewportBuffer.fbo);

	glViewport(0, 0, viewPortDimensions.width, viewPortDimensions.height);

	glClearColor(viewportBGColour.x, viewportBGColour.y, viewportBGColour.z, viewportBGColour.w);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//test draw supertriangle
	if (nodes.size() > 2) 
	{
		glUseProgram(testSuperTriShader.program);
		glBindVertexArray(layers[10].meshData.vertexArrayObject);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, layers[10].meshData.vertexArrayElementObject);
		glDrawArrays(GL_TRIANGLES, 0, 3);
	}

	//The nodes and triangles share the same vertices, we use the same vertex object for both, but we glDrawArrays the point, and for
	//the elements we use glDrawElements.
	glBindVertexArray(layers[1].meshData.vertexArrayObject);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, layers[1].meshData.vertexArrayElementObject);
	
	//draw triangles
	glUseProgram(triangleShader.program);
	glDrawElements(GL_TRIANGLES, renderTrisIndicesCount, GL_UNSIGNED_INT, 0);
	
	//draw nodes
	glUseProgram(pointShader.program);
	glDrawArrays(GL_POINTS, 0, renderVertsCount);

	//restore default frame buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

//IMGUI drawing commands, passes the texture of the offscreen buffer to IMGUI to be painted in the viewport window.
void DrawViewport() //IMGUI painting commands for viewport window
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	ImGui::SetNextWindowPos(ImVec2(viewPortDimensions.positionX, viewPortDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(viewPortDimensions.width, viewPortDimensions.height), ImGuiCond_Always);

	ImGui::Begin("Viewport", NULL, windowFlags);

	ImVec2 pos = ImGui::GetCursorScreenPos();

	ImGui::GetWindowDrawList()->AddImage(
		(void*)viewportBuffer.texture, //texture to render (from an offscreen buffer)
		ImVec2(viewPortDimensions.positionX, viewPortDimensions.positionY), //NW Corner
		ImVec2(viewPortDimensions.positionX + viewPortDimensions.width, viewPortDimensions.positionY + viewPortDimensions.height), //SE corner
		ImVec2(0, 1), ImVec2(1, 0)); //UVs

	ImGui::End();
}

bool InitViewport()
{
	viewPortDimensions.width = minViewportWidth;
	viewPortDimensions.height = minViewportHeight;

	LogMan::Log("Initializaing viewport");
	
	if (!SetupOffScreenBuffer(&viewportBuffer, viewPortDimensions.width, viewPortDimensions.height))
		return FAILED_VIEWPORT_CREATE;

	layers.insert({ 1, Layer("MainMesh", MeshData(), 1) });
	//layers.insert({ 2, Layer("Triangles", MeshData(), 2) });
	//layers.insert({ 5, Layer("Test", MeshData(), 5) });
	
	for (auto it = layers.begin(); it != layers.end(); ++it)
		SetupMesh(&(it->second.meshData), NULL, 0, NULL, 0); //startup empty

	if ( !SetupShaders(&pointShader, CircleShader::vertex_shader_text, CircleShader::fragment_shader_text) ||
		!SetupShaders(&triangleShader, TestTriangle::vertex_shader_text, TestTriangle::fragment_shader_text) )
	{
		LogMan::Log("Error compiling shaders!", LOG_ERROR);
		//shader compilation failure is not [very] critical. Most GPUs offer fallback shaders in this case anyway,
		//although they may not show much, but at least the program as a whole won't break. No need to terminate
		//the entire program just for this. Although more graceful handling of the issue should be implemented.
		//return false;
	}
	
	//Check and push any errors that may have resulted thus far to the log.
	GLErrorCheck();

	TestSetupSuperTriangleRender(); //test
	return true;
}

//To be called when the dimensions of the viewport change (e.g. when windows are resized).
void UpdateViewport()
{
	UpdateOffScreenBuffer(&viewportBuffer, viewPortDimensions.width, viewPortDimensions.height);
	UpdateViewBounds();
	UpdateCoordinateProjectionParameters();
	UpdateContent();
}

//Recomputes the viewport-space positions of objects to be drawn.
void UpdateContent()
{
	UpdateNodes();
	UpdateTriangles();

	UpdateMesh(&(layers[1].meshData), nodesMeshVerts, renderVertsCount, trianglesIndices, renderTrisIndicesCount);
	TestUpdateSuperTriangle();
}

//Recomputes the viewport-space positions of nodes as stored in the ModelInterface. If nodes count is less than 1, this function \
only clears the viewport mesh container of the nodes.
void UpdateNodes()
{
	CLEAR_ARRAY(nodesMeshVerts)
	renderVertsCount = nodes.size() * 3;

	if (nodes.size() < 1)
		return;

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

//Fills the triangleIndices array with the vertexIDs from each triangle object in triangles map stored in ModelInterface. \
If triangles count is less than 1, only clears the indices array. 
void UpdateTriangles()
{
	CLEAR_ARRAY(trianglesIndices)
	renderTrisIndicesCount = triangles.size() * 3;

	if (triangles.size() < 1)
		return;

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


//To be called whenever viewport is resized. Recomputes the NE corner based on the delta of the new position (assumed to be updated in \
viewPortDimensions) and the previous viewport (Assumed stored in lastViewportSize). SW corner does not change. Recomputes aspectRatio.
void UpdateViewBounds()
{
	//scale doesn't change when simply resizing viewport. (i.e. resizing would expose more area to paint on, but won't change size of elements.
	aspectRatio = static_cast<double>(viewPortDimensions.width) / static_cast<double>(viewPortDimensions.height);
	
	Vector2 changeInView = viewPortDimensions.Dimension() - lastViewportSize;	
	viewBounds[1].x += changeInView.x * scale;
	viewBounds[1].y += changeInView.y * scale / aspectRatio;
}

//Forces the viewportbounds to specific coordinates. SW bound is set as provided, NW corner is clamped to be at least \
MIN_VIEWPORT_DELTA north and east (each) of the SW corner. New NW corner is adjusted to maintain current window aspect ratio. \
Updates scale and delta.
void SetViewBounds(Vector2 swCorner, Vector2 nwCorner)
{
	viewBounds[0] = swCorner;
	viewBounds[1] = Vector2(Max(nwCorner.x, swCorner.x + MIN_VIEWPORT_DELTA), Max(nwCorner.y, swCorner.y + MIN_VIEWPORT_DELTA));

	//adjust NW bound to maintain the current aspect ratio.
	delta = viewBounds[1] - viewBounds[0]; //unadjusted delta
	delta = Vector2(Max(delta.x, delta.y * aspectRatio),
					Max(delta.y, delta.x / aspectRatio));

	viewBounds[1] = viewBounds[0] + delta;

	scale = delta.x / static_cast<double>(viewPortDimensions.width);
}

void UpdateCoordinateProjectionParameters()
{
	delta = viewBounds[1] - viewBounds[0];
	aspectRatio = static_cast<double>(viewPortDimensions.width) / static_cast<double>(viewPortDimensions.height);
	scale = delta.x / static_cast<double>(viewPortDimensions.width);
}

//returns normalized coordinates (-1.0f to 1.0f) for a supplied point depending on the current projection parameters 
Vector2 NormalizeCoordinates(Vector2 & point)
{
	/*Vector2 normPos((2.0f - 2.0f * margin) * (point.x - viewBounds[0].x) / (delta.x) - 1.0f + margin,
					(2.0f - 2.0f * margin) * (point.y - viewBounds[06].y) / (delta.y) - 1.0f + margin);*/
	/*Vector2 normPos(2.0f * (point.x - viewBounds[0].x) / (delta.x) - 1.0f,
					2.0f * (point.y - viewBounds[0].y) / (delta.y) - 1.0f);*/
	Vector2 normPos(2.0f * ((point.x - viewBounds[0].x) / scale) / viewPortDimensions.width  - 1.0f,
					2.0f * ((point.y - viewBounds[0].y) / scale) / viewPortDimensions.height - 1.0f);
	
	normPos = normPos;
	return normPos;
}

//Shifts the viewBounds by posDelta. Calls UpdateContent before returning.
void PanView(Vector2 posDelta)
{
	SetViewBounds(viewBounds[0] + posDelta, viewBounds[1] + posDelta);
	UpdateContent();
}

//Stuff for debug display of the supertriangle.
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