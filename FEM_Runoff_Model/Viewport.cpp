#include "Viewport.hpp"

OffScreenBuffer viewportBuffer;

Vector2D delta;
double scale = 0.5f;
double screenAspectRatio = 1.0f;
float scaleChangeTicks = 0.01f;
float * nodesMeshVerts = NULL;
float * boundaryLineVerts = NULL;
unsigned int * trianglesIndices = NULL;
int renderVertsCount = 0;
int renderTrisIndicesCount = 0;
int renderBoundVertsCount = 0;

std::unordered_map <int, Layer>  layers;
Shader triangleShader, pointShader, lineShader;

Vector2D lastViewportSize;
Vector2D viewBounds[2];
Vector2D currenViewportHoverPos;
Vector2D currenViewportHoverPosPixels; //in local pixel screenspace of viewport
bool isHoveringViewport = false;

bool showElementView = false;
Triangle const * selectedElement = NULL;
//Cached strings of element's details, to avoid creating strings every loop.
std::string elementViewName = "";
std::string elementViewVerts = "";
std::string elementViewArea = "";

ToolMode activeTool = ToolMode::None;
bool isGrabbingNode = false;
size_t activeNode = 0;

float currentPointSize = 10.0f; //TODO the shader uses "float" for point size, but the gl_PointSize is actually in pixels, which should be integer. Check this.

#pragma region Forward declarations
//Several (if not most) don't really need to be declared here, but it was easier to copy and paste them from the header (where they were initially).
void SetupMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount);
void UpdateMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount);
bool SetupShaders(Shader * targetGLData, char const * vertexShaderSource, char const * fragmentShaderSource);
bool SetupOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY);
void UpdateOffScreenBuffer(OffScreenBuffer * buffer, int sizeX, int sizeY);
void RenderViewport();
void UpdateContent();
void UpdateNodes();
void UpdateTriangles();
void UpdateWatershed();
void UpdateViewBounds();
void UpdateCoordinateProjectionParameters();
Vector2D NormalizeCoordinates(Vector2D const & point);
Vector2D ScreenToWorldSpace(Vector2D const & screenCoordinates);
void PanView(Vector2D posDelta);
#pragma endregion

//Create vertex buffer (VBO), vertex array object (VAO) and vertex array element object, store id in targetGLData struct. \
if verticesCount is greater than or equal to 3, UpdateMesh() is called.
void SetupMesh(MeshData * targetGLData, float const * mesh, unsigned int verticesCount, unsigned int const * indices, unsigned int indexCount)
{
	glGenBuffers(1, &(targetGLData->vertexBufferObject));
	glGenVertexArrays(1, &(targetGLData->vertexArrayObject));
	glGenBuffers(1, &(targetGLData->vertexArrayElementObject));
	
	if (verticesCount >= 1)
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
	glBufferData(GL_ARRAY_BUFFER, verticesCount * 3 * sizeof(float), mesh, GL_STATIC_DRAW);
	
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

void SetActiveShaderDiffuse(Colour const & colour)
{
	GLuint diffuseColLoc = glGetUniformLocation(triangleShader.program, "diffuseCol");
	glUniform4fv(diffuseColLoc, 1, colour.Array());
}

//Note: Switches active program to pointShader.program. Does not reset to existing program before calling.
void SetPointSize(float size)
{
	glUseProgram(pointShader.program);
	GLuint pointSizeLoc = glGetUniformLocation(pointShader.program, "pointSize");
	glUniform1f(pointSizeLoc, size);
	currentPointSize = size;
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

ToolMode GetActiveTool()
{
	return activeTool;
}

void SwitchActiveTool(ToolMode newTool)
{
	isGrabbingNode = false;
	activeTool = newTool;
}

void HandleMousePan()
{
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 * posAtClick = io.MouseClickedPos;
	if (ImGui::IsMouseDragging(2) && viewportDimensions.Contains(posAtClick[2]))
	{
		ImVec2 mouseDelta = io.MouseDelta;
		PanView(Vector2D(-1.0f * mouseDelta.x, mouseDelta.y) * scale);
	}
}

void HandleMouseZoom()
{
	ImGuiIO& io = ImGui::GetIO();
	float scroll = io.MouseWheel * -1.0f;
	if (scroll != 0.0f)
	{
		//Compute a temporary scale and bounds to test.
		double newScale = scale + scroll * scaleChangeTicks;
		Vector2D newBounds[2] = {	Vector2D(	currenViewportHoverPos.x - currenViewportHoverPosPixels.x * newScale,
												currenViewportHoverPos.y - (viewportDimensions.height - currenViewportHoverPosPixels.y) * newScale),
									Vector2D(	currenViewportHoverPos.x + (viewportDimensions.width - currenViewportHoverPosPixels.x) * newScale,
												currenViewportHoverPos.y + currenViewportHoverPosPixels.y * newScale) };

		//Test is that if span in either axes of the new bounds is less than MIN_VIEWPORT_DELTA, we don't zoom.
		Vector2D diff = newBounds[1] - newBounds[0];
		if (diff.x < MIN_VIEWPORT_DELTA || diff.y < MIN_VIEWPORT_DELTA)
			return;

		viewBounds[0] = newBounds[0];
		viewBounds[1] = newBounds[1];
		scale = newScale;

		UpdateCoordinateProjectionParameters();
		UpdateContent();
	}
}

bool SelectHoveredElement()
{
	selectedElement = GetElementContainingPoint(currenViewportHoverPos);
	
	if (selectedElement != NULL)
	{
		elementViewName = std::move(std::string("Element: " + std::to_string(selectedElement->id)).c_str());
		elementViewVerts = std::move(std::string("Vertices: " + std::to_string(selectedElement->VertexID(0)) + ", "
			+ std::to_string(selectedElement->VertexID(1)) + ", "
			+ std::to_string(selectedElement->VertexID(2))).c_str());
		elementViewArea = std::move(std::string("Area: " + std::to_string(selectedElement->Area())).c_str());
		return true;	
	}
	else
	{
		elementViewName = "";
		elementViewVerts = "";
		elementViewArea = "";
		return false;
	}
}

void HandleMouseLeftClick()
{
	//check that we're clicking on
	if (ImGui::IsMouseClicked(0))
	{
		switch (activeTool)
		{
			case None:
			{
				return;
			}
			case ViewNode:
			{
				LogMan::Log("Viewing nodes not implemented yet.", LOG_WARN);
				return;
			}
			case ViewElement:
			{
				if (SelectHoveredElement())
				{
					ImGui::OpenPopup("elementView");
					showElementView = true;
				}
				else
					showElementView = false;
				return;
			}
			case MoveNode:
			{
				if (isGrabbingNode)
				{
					//place node
					//LogMan::Log("Placing node: " + std::to_string(activeNode) + " at " + std::to_string(currenViewportHoverPos.x) + ", " + std::to_string(currenViewportHoverPos.y)); //test
					UpdateNode(activeNode, currenViewportHoverPos);
					isGrabbingNode = false;
				}
				else
				{
					for (int i = 0; i < GetNodes().size(); i++)
					{
						//check of within range of node
						//our range is a circle of centre = node's position and radius = currentPointSize.
						if (currenViewportHoverPos.WithinCircle(GetNodes()[i], currentPointSize * scale))
						{
							//grab node
							//LogMan::Log("Grabbed node: " + std::to_string(activeNode)); //test
							isGrabbingNode = true;
							activeNode = i;
							break;
						}
					}
				}
				return;
			}
			default:
			{
				return;
			}
		}
	}
}

//MUST BE CALLED BEFORE Viewport input handling. \
Only updates when position when mouse is over viewport, otherwise isHoveringViewport set to false.
void UpdateMouseHoverPosition()
{
	ImGuiIO& io = ImGui::GetIO();
	ImVec2 mousePos = io.MousePos;
	
	if (isHoveringViewport = viewportDimensions.Contains(mousePos))
	{
		currenViewportHoverPosPixels = viewportDimensions.LocalPosFromGlobal(mousePos);
		currenViewportHoverPos = ScreenToWorldSpace(currenViewportHoverPosPixels);
	}
	else
		isHoveringViewport = false;
}

void UpdateGrabbedNode()
{
	if (!isGrabbingNode)
		return;
	
	//recompute the screenspace position of the grabbed node.
	//Note that this won't affect the original position stored in nodes vector...
	Vector2D relativePos = NormalizeCoordinates(currenViewportHoverPos);

	nodesMeshVerts[activeNode * 3] = relativePos.x;
	nodesMeshVerts[(activeNode * 3) + 1] = relativePos.y;
	nodesMeshVerts[(activeNode * 3) + 2] = 0.0f;
	UpdateMesh(&(layers[1].meshData), nodesMeshVerts, renderVertsCount, trianglesIndices, renderTrisIndicesCount);
}

//called by RenderViewport Only.
void RenderMesh()
{
	//Switch to FEM mesh
		//The nodes and triangles share the same vertices, we use the same vertex object for both, but we glDrawArrays the point, and for
		//the elements we use glDrawElements.
	glBindVertexArray(layers[1].meshData.vertexArrayObject);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, layers[1].meshData.vertexArrayElementObject);

	//draw triangles
	glUseProgram(triangleShader.program);
	SetActiveShaderDiffuse(COLOUR_RED);
	glDrawElements(GL_TRIANGLES, renderTrisIndicesCount, GL_UNSIGNED_INT, 0);
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //Wireframe only
	SetActiveShaderDiffuse(COLOUR_BLACK);
	glDrawElements(GL_TRIANGLES, renderTrisIndicesCount, GL_UNSIGNED_INT, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); //return to normal filled render mode.

	//draw nodes
	glUseProgram(pointShader.program);
	SetActiveShaderDiffuse(COLOUR_BLACK);
	glDrawArrays(GL_POINTS, 0, renderVertsCount);
}

//Renders the viewport content to the offscreen buffer.
void RenderViewport() //Renders viewport content to an offscreen buffer.
{
	//Objects in the back are rendered first.

	glBindFramebuffer(GL_FRAMEBUFFER, viewportBuffer.fbo);

	glViewport(0, 0, viewportDimensions.width, viewportDimensions.height);

	glClearColor(viewportBGColour.x, viewportBGColour.y, viewportBGColour.z, viewportBGColour.w);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	RenderMesh();

	//Draw watershed boundary
	glBindVertexArray(layers[2].meshData.vertexArrayObject);
	glUseProgram(triangleShader.program);
	SetActiveShaderDiffuse(COLOUR_MAGENTA);
	glDrawArrays(GL_LINE_LOOP, 0, renderBoundVertsCount);

	//restore default frame buffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

//Renders viewport then draws to screen IMGUI drawing commands, passes the texture of the offscreen buffer to IMGUI to \
be painted in the viewport window.
void DrawViewport()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;

	ImGui::Begin("Viewport", NULL, windowFlags);
	
	UpdateMouseHoverPosition();
	if (isHoveringViewport)
	{
		HandleMousePan();
		HandleMouseZoom();
		HandleMouseLeftClick();
	}
	UpdateGrabbedNode();

	RenderViewport();

	ImVec2 pos = ImGui::GetCursorScreenPos();
	
	ImGui::SetWindowPos(ImVec2(viewportDimensions.positionX, viewportDimensions.positionY), ImGuiCond_Always);
	ImGui::SetWindowSize(ImVec2(viewportDimensions.width, viewportDimensions.height), ImGuiCond_Always);

	ImGui::GetWindowDrawList()->AddImage(
		(void*)viewportBuffer.texture, //texture to render (from an offscreen buffer)
		ImVec2(viewportDimensions.positionX, viewportDimensions.positionY), //NW Corner
		ImVec2(viewportDimensions.positionX + viewportDimensions.width, viewportDimensions.positionY + viewportDimensions.height), //SE corner
		ImVec2(0, 1), ImVec2(1, 0)); //UVs

	if (showElementView) //redundant?
	{
		if (ImGui::BeginPopup("elementView"))
		{
			ImGui::Text(elementViewName.c_str());
			ImGui::Text(elementViewVerts.c_str());
			ImGui::Text(elementViewArea.c_str());
			ImGui::EndPopup();
		}
	}

	ImGui::End();
}

bool InitViewport()
{
	viewportDimensions.width = 1;
	viewportDimensions.height = 1;

	LogMan::Log("Initializaing viewport");
	
	if (!SetupOffScreenBuffer(&viewportBuffer, viewportDimensions.width, viewportDimensions.height))
		return FAILED_VIEWPORT_CREATE;

	layers.insert({ 1, Layer("MainMesh", MeshData(), 1) });
	layers.insert({ 2, Layer("WatershedBoundary", MeshData(), 2) });
	//layers.insert({ 2, Layer("Triangles", MeshData(), 2) });
	//layers.insert({ 5, Layer("Test", MeshData(), 5) });

	for (auto it = layers.begin(); it != layers.end(); ++it)
		SetupMesh(&(it->second.meshData), NULL, 0, NULL, 0); //startup empty

	if ( !SetupShaders(&pointShader, PointShader::vertex_shader_text, PointShader::fragment_shader_text) ||
		!SetupShaders(&triangleShader, PolygonShader::vertex_shader_text, PolygonShader::fragment_shader_text) )
	{
		LogMan::Log("Error compiling shaders!", LOG_ERROR);
		//shader compilation failure is not [very] critical. Most GPUs offer fallback shaders in this case anyway,
		//although they may not show much, but at least the program as a whole won't break. No need to terminate
		//the entire program just for this. Although more graceful handling of the issue should be implemented.
		//return false;
	}

	//Set initial point render size for point shader
	SetPointSize(currentPointSize);

	//Check and push any errors that may have resulted thus far to the log.
	GLErrorCheck();

	return true;
}

//To be called when the dimensions of the viewport change (e.g. when windows are resized).
void UpdateViewport()
{
	UpdateOffScreenBuffer(&viewportBuffer, viewportDimensions.width, viewportDimensions.height);
	UpdateViewBounds();
	//UpdateCoordinateProjectionParameters();
	UpdateContent();
}

//Recomputes the viewport-space positions of objects to be drawn.
void UpdateContent()
{
	UpdateNodes();
	UpdateTriangles();
	UpdateWatershed();

	UpdateMesh(&(layers[1].meshData), nodesMeshVerts, renderVertsCount, trianglesIndices, renderTrisIndicesCount);
	UpdateMesh(&(layers[2].meshData), boundaryLineVerts, renderBoundVertsCount, NULL, 0);
}

//Recomputes the viewport-space positions of nodes as stored in the ModelInterface. If nodes count is less than 1, this function \
only clears the viewport mesh container of the nodes.
void UpdateNodes()
{
	//TODO this function is identical to UpdatWatershed, and probably all similar future functions. Unify in a one method that takes\
	the variable params, and call for each set of params from UpdateContent.

	CLEAR_ARRAY(nodesMeshVerts);
	
	auto nodes = GetNodes();

	renderVertsCount = nodes.size();

	if (renderVertsCount < 1)
		return;
	
	nodesMeshVerts = new float[renderVertsCount * 3];

	int counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		Vector2D relativePos = NormalizeCoordinates(*it);
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
	CLEAR_ARRAY(trianglesIndices);
	
	auto triangles = GetTriangles();
	renderTrisIndicesCount = triangles.size() * 3;
	
	if (triangles.size() < 1)
		return;

	trianglesIndices = new unsigned int[renderTrisIndicesCount];

	int counter = 0;
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		trianglesIndices[counter] = it->second.VertexID(0);
		trianglesIndices[counter + 1] = it->second.VertexID(1);
		trianglesIndices[counter + 2] = it->second.VertexID(2);
		counter += 3;
	}
}

void UpdateWatershed()
{
	CLEAR_ARRAY(boundaryLineVerts);
	
	auto shedBoundary = GetWatershedBoundary();
	renderBoundVertsCount = shedBoundary.size();
	if (renderBoundVertsCount < 2)
	{
		renderBoundVertsCount = 0;
		return;
	}

	boundaryLineVerts = new float[renderBoundVertsCount * 3];
	int counter = 0;

	for (auto it = shedBoundary.begin(); it != shedBoundary.end(); ++it)
	{
		Vector2D relativePos = NormalizeCoordinates(*it);
		boundaryLineVerts[counter] = relativePos.x;
		boundaryLineVerts[counter + 1] = relativePos.y;
		boundaryLineVerts[counter + 2] = 0.0f;
		counter += 3;
	}
}

//To be called whenever viewport is resized. Recomputes the NE corner based on the delta of the new position (assumed to be updated in \
viewportDimensions) and the previous viewport (Assumed stored in lastViewportSize). SW corner does not change. Recomputes screenAspectRatio. \
Imples UpdateCoordinateProjectionParameters()
void UpdateViewBounds()
{
	//scale doesn't change when simply resizing viewport. (i.e. resizing would expose more area to paint on, but won't change size of elements.
	screenAspectRatio = static_cast<double>(viewportDimensions.width) / static_cast<double>(viewportDimensions.height);
	
	Vector2D changeInView = static_cast<Vector2D>(viewportDimensions.Dimension()) - lastViewportSize;	
	viewBounds[1].x += changeInView.x * scale;
	viewBounds[1].y += changeInView.y * scale / screenAspectRatio;

	UpdateCoordinateProjectionParameters();
}

//Forces the viewportbounds to specific coordinates. SW bound is set as provided, NW corner is clamped to be at least \
MIN_VIEWPORT_DELTA north and east (each) of the SW corner. New NW corner is adjusted to maintain current window aspect ratio. \
Updates scale and delta.
void SetViewBounds(Vector2D const & swCorner, Vector2D const & neCorner)
{
	//Sanitize the viewbounds before doing anything (viewBounds can't be flipped, and can't be the same).
	viewBounds[0] = swCorner;
	viewBounds[1] = Vector2D(Max(static_cast<double>(neCorner.x), swCorner.x + MIN_VIEWPORT_DELTA),
							Max(static_cast<double>(neCorner.y), swCorner.y + MIN_VIEWPORT_DELTA));

	//adjust NE bound to maintain the current aspect ratio.
	delta = viewBounds[1] - viewBounds[0]; //unadjusted delta
	delta = Vector2D(Max(static_cast<double>(delta.x), delta.y * screenAspectRatio),
					Max(static_cast<double>(delta.y), delta.x / screenAspectRatio));

	viewBounds[1] = viewBounds[0] + delta;
	UpdateCoordinateProjectionParameters(); //Although it reduntatly recomputes delta...
	UpdateContent();
}

void CentreOnObject(Vector2D const & swCorner, Vector2D const & neCorner)
{
	//Sanitize the viewbounds before doing anything (viewBounds can't be flipped, and can't be the same).
	viewBounds[0] = swCorner;
	viewBounds[1] = Vector2D(Max(static_cast<double>(neCorner.x), swCorner.x + MIN_VIEWPORT_DELTA),
							Max(static_cast<double>(neCorner.y), swCorner.y + MIN_VIEWPORT_DELTA));

	//Code bellow will handle fixing aspect ratio.	
	viewBounds[0].x -= VIEW_CENTRE_PADDING * scale;
	viewBounds[0].y -= VIEW_CENTRE_PADDING * scale;
	viewBounds[1].x += VIEW_CENTRE_PADDING * scale;
	viewBounds[1].y += VIEW_CENTRE_PADDING * scale;


	//adjust NE bound to maintain the current aspect ratio.
	delta = viewBounds[1] - viewBounds[0]; //unadjusted delta
	delta = Vector2D(Max(static_cast<double>(delta.x), delta.y * screenAspectRatio),
					Max(static_cast<double>(delta.y), delta.x / screenAspectRatio));

	viewBounds[1] = viewBounds[0] + delta;

	Vector2D addedPadding = viewBounds[1] - neCorner;

	viewBounds[0] = viewBounds[0] - addedPadding * 0.5;
	viewBounds[1] = viewBounds[1] - addedPadding * 0.5;

	UpdateCoordinateProjectionParameters(); //Although it reduntatly recomputes delta...
	UpdateContent();
}

void UpdateCoordinateProjectionParameters()
{
	delta = viewBounds[1] - viewBounds[0];
	screenAspectRatio = static_cast<double>(viewportDimensions.width) / static_cast<double>(viewportDimensions.height);
	scale = delta.x / static_cast<double>(viewportDimensions.width);
	scaleChangeTicks = 0.1 * scale;
}

//returns normalized coordinates (-1.0f to 1.0f) for a supplied point depending on the current projection parameters 
Vector2D NormalizeCoordinates(Vector2D const & point)
{
	Vector2D normPos(2.0f * ((point.x - viewBounds[0].x) / scale) / viewportDimensions.width  - 1.0f,
					2.0f * ((point.y - viewBounds[0].y) / scale) / viewportDimensions.height - 1.0f);
	
	return normPos;
}

Vector2D ScreenToWorldSpace(Vector2D const & screenCoordinates)
{
	return Vector2D(viewBounds[0].x + (screenCoordinates.x * scale),
					viewBounds[1].y - (screenCoordinates.y * scale));
}

//Shifts the viewBounds by posDelta. Calls UpdateContent before returning.
void PanView(Vector2D posDelta)
{
	SetViewBounds(viewBounds[0] + posDelta, viewBounds[1] + posDelta);
	UpdateContent();
}