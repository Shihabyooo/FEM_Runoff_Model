#include "GUI_Handler.hpp"
#include "GUI_Requirements.hpp"

const int minMainWinWidth = 800, minMainWinHeight = 600;
const int minViewportWidth = 800, minViewportHeight = 600;

int mainWinWidth = 1280, mainWinHeight = 720;
int viewportWidth = 800, viewportHeight = 600;

GLFWwindow * mainWindow;
GLData viewportGLData;
OffScreenBuffer viewportBuffer;
std::list<OffScreenBuffer> offScreenBuffersList;

ImVec4 mainBGColour = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
ImVec4 viewportBGColour = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

WindowDimensions leftPaneDimensions, logPaneDimensions, viewPortDimensions;

void OnGLFWError(int error, const char* description) //Callback
{
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

void OnMainWindowSizeChange(GLFWwindow * window, int width, int height)
{
	mainWinWidth = width;
	mainWinHeight = height;

	RecomputeWindowElementsDimensions();// width, height);
}

void InitializeDearIMGUI(const char * glslVersion) //code copied directly from official examples
{
	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsClassic();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(mainWindow, true);
	ImGui_ImplOpenGL3_Init(glslVersion);

	ImVec4 mainBGColour = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
}

bool InitializeMainWindow() //implies InitializeDearIMGUI()
{
	//set error callback
	glfwSetErrorCallback(OnGLFWError);

	if (!glfwInit())
	{
		std::cout << "Failed to initialize GLFW!" << std::endl;
		return false;
	}

	//attempt to use GL 3.0
	const char* glslVersion = "#version 130"; //needed by imgui ogl impl
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	glfwWindowHint(GLFW_SAMPLES, 4); //use 4 samples for AA

	mainWindow = glfwCreateWindow(mainWinWidth, mainWinHeight, PROGRAM_NAME, NULL, NULL);
	if (mainWindow == NULL)
	{
		std::cout << "Failed to create main window!" << std::endl;
		glfwTerminate();
		return false;
	}

	glfwMakeContextCurrent(mainWindow);
	glfwSwapInterval(1);
	glfwSetWindowSizeCallback(mainWindow, OnMainWindowSizeChange);
	glfwSetWindowSizeLimits(mainWindow, minMainWinWidth, minMainWinHeight, GLFW_DONT_CARE, GLFW_DONT_CARE); //set min screen dims
	//glewExperimental = GL_TRUE;
	glewInit();

	InitializeDearIMGUI(glslVersion);

	return true;
}

void TerminateMainWindow()
{
	glfwDestroyWindow(mainWindow);
	glfwTerminate();
}

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

void RecomputeWindowElementsDimensions()// int newMainWinWidth, int newMainWinHeight)
{
	//Set fixed width for left pane.
	//Set fixed height for log pane.
	//left pane height = mainWindow height.
	//log pane width = mainwindow width - fixed left pane width.
	//viewport takes remainin area.
	int fixedLeftPaneWidth = 250;
	int fixedLogPaneHeight = 100;

	leftPaneDimensions = WindowDimensions(0, 0, fixedLeftPaneWidth, mainWinHeight);
	logPaneDimensions = WindowDimensions(leftPaneDimensions.width, mainWinHeight - fixedLogPaneHeight, mainWinWidth - fixedLeftPaneWidth, fixedLogPaneHeight);
	viewportWidth = mainWinWidth - fixedLeftPaneWidth;
	viewportHeight = mainWinHeight - fixedLogPaneHeight;
	viewPortDimensions = WindowDimensions(leftPaneDimensions.width, 0, viewportWidth, viewportHeight);

	//update buffer for viewport to use new dimensions.
	UpdateOffScreenBuffer(&viewportBuffer, viewportWidth, viewportHeight);
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

//TODO move window/viewport rendering/painting function to their own header/source files.
char gridNodes[260] = "Path to grid nodes file.";
char demFilePath[260] = "Path to grid nodes file.";

void DrawLeftPane()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	ImGui::SetNextWindowPos(ImVec2(leftPaneDimensions.positionX, leftPaneDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(leftPaneDimensions.width, leftPaneDimensions.height), ImGuiCond_Always);

	ImGui::Begin("Input", NULL, windowFlags);
	ImGui::PushItemWidth(ImGui::GetFontSize() * -12); //Use fixed width for labels (by passing a negative value), the rest goes to widgets. We choose a width proportional to our font size.
	
	if (ImGui::BeginMenuBar())
	{
		if (ImGui::BeginMenu("File"))
			//ImGui::MenuItem("New", NULL, &newFile);
			ImGui::EndMenu();
		if (ImGui::BeginMenu("Edit"))
			ImGui::EndMenu();
		if (ImGui::BeginMenu("Help"))
			ImGui::EndMenu();
		ImGui::EndMenuBar();
	}
	
	ImGui::Text(PROGRAM_NAME);
	ImGui::Separator();
	ImGui::NewLine();

	//Geometry data
	ImGui::Text("Mesh Nodes");
	ImGui::InputText("Mesh Nodes", gridNodes, IM_ARRAYSIZE(gridNodes));

	if (ImGui::Button("Browse for geometry directory"))
	{
		//TODO spawn file browser here
	}

	/*DrawFileBrowser();
	DrawFileList(geometryFilePath, &geometryNames, &selectedGeometry, DataType::geometry, false, GEOMETRY_LIST_ID);
	ImGui::NewLine();*/

	//DEM
	ImGui::Text("DEM");
	ImGui::InputText("DEM File Path", demFilePath, IM_ARRAYSIZE(demFilePath));

	if (ImGui::Button("Browse for DEM directory"))
	{

	}
	////DrawFileBrowser(); //The call is already made above...
	//DrawFileList(demFilePath, &demNames, &selectedDEM, DataType::dem, true, DEM_LIST_ID);
	//ImGui::NewLine();
		
	//Other input

	ImGui::NewLine();
	ImGui::NewLine();
	if (ImGui::Button("Run Simulation!", ImVec2(100, 50)))
	{

	}

	ImGui::End();
}

void DrawLogPane()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	ImGui::SetNextWindowPos(ImVec2(logPaneDimensions.positionX, logPaneDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(logPaneDimensions.width, logPaneDimensions.height), ImGuiCond_Always);
	
	ImGui::Begin("Log", NULL, windowFlags);

	ImGui::PushItemWidth(ImGui::GetFontSize() * -12); //Use fixed width for labels (by passing a negative value), the rest goes to widgets. We choose a width proportional to our font s

	ImGui::TextWrapped("Log goes here");
	ImGui::End();
}

void RenderViewport() //Renders viewport content to an offscreen buffer.
{
	//Draw to offscreen buffer.
	glBindFramebuffer(GL_FRAMEBUFFER, viewportBuffer.fbo);
	int minDimension = viewportWidth > viewportHeight ? viewportHeight : viewportWidth;
	glViewport(0, 0, minDimension, minDimension);
	glClearColor(viewportBGColour.x, viewportBGColour.y, viewportBGColour.z, viewportBGColour.w);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_DEPTH_TEST);

	glUseProgram(viewportGLData.program);
	glBindVertexArray(viewportGLData.vertexArrayObject);
	glDrawArrays(GL_TRIANGLES, 0, 3);

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

int MainUILoop()
{
	bool showDemo = true;
	while (!glfwWindowShouldClose(mainWindow))
	{
		glfwPollEvents();

		glViewport(0, 0, mainWinWidth, mainWinHeight); //update viewport to current window size

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		glClearColor(mainBGColour.x * mainBGColour.w, mainBGColour.y * mainBGColour.w, mainBGColour.z * mainBGColour.w, mainBGColour.w);
		glClear(GL_COLOR_BUFFER_BIT);

		ImGui::ShowDemoWindow(&showDemo);

		DrawLeftPane();
		DrawLogPane();

		RenderViewport();
		DrawViewport();

		ImGui::Render();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		glfwSwapBuffers(mainWindow);
	}
	
	DeleteAllOffscreenBuffers();
	TerminateMainWindow();

	return SUCCESS;
}

int StartUI()
{
	if (!InitializeMainWindow())
	{
		std::cout << "Failed to initialize main window" << std::endl;
		return FAILED_MAIN_WINDOW_INITIALIZE;
	}

	SetupMesh(&viewportGLData, TestTriangle::vertices, sizeof(TestTriangle::vertices));
	SetupShaders(&viewportGLData, const_cast<char *>(TestTriangle::vertex_shader_text), const_cast<char *>(TestTriangle::fragment_shader_text));

	if (!SetupOffScreenBuffer(&viewportBuffer, viewportWidth, viewportHeight))
	{
		TerminateMainWindow();
		return FAILED_VIEWPORT_CREATE;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0); //restore default framebuffer

	RecomputeWindowElementsDimensions();// mainWinWidth, mainWinHeight);

	return MainUILoop();
}