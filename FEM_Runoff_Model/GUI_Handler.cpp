#pragma once
#include "GUI_Handler.hpp"
#include "GUI_Requirements.hpp"
#include "Viewport.hpp"
#include "LogManager.hpp"

int mainWinWidth = 1280, mainWinHeight = 720;
int viewportWidth = 800, viewportHeight = 600;

GLFWwindow * mainWindow;

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
	const char* glslVersion = "#version 330"; //needed by imgui ogl impl
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
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

void RecomputeWindowElementsDimensions()// int newMainWinWidth, int newMainWinHeight)
{
	//Set fixed width for left pane.
	//Set fixed height for log pane.
	//left pane height = mainWindow height.
	//log pane width = mainwindow width - fixed left pane width.
	//viewport takes remainin area.
	int fixedLeftPaneWidth = 250;
	int fixedLogPaneHeight = 200;

	leftPaneDimensions = WindowDimensions(0, 0, fixedLeftPaneWidth, mainWinHeight);
	logPaneDimensions = WindowDimensions(leftPaneDimensions.width, mainWinHeight - fixedLogPaneHeight, mainWinWidth - fixedLeftPaneWidth, fixedLogPaneHeight);
	viewportWidth = mainWinWidth - fixedLeftPaneWidth;
	viewportHeight = mainWinHeight - fixedLogPaneHeight;
	viewPortDimensions = WindowDimensions(leftPaneDimensions.width, 0, viewportWidth, viewportHeight);

	//update buffer for viewport to use new dimensions.
	//UpdateOffScreenBuffer(&viewportBuffer, viewportWidth, viewportHeight);
	UpdateViewport();
}

//TODO move window/viewport rendering/painting function to their own header/source files.
//char meshNodes[260] = "Path to grid nodes file.";
char meshNodes[260] = "Test_Mesh_Nodes_Grid.csv";
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
	ImGui::PushItemWidth(-1);
	ImGui::InputText("Mesh Nodes", meshNodes, IM_ARRAYSIZE(meshNodes));
	ImGui::PopItemWidth();
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
		std::string nodePath(meshNodes);
		TestSimulate(nodePath);
		UpdateViewport();
	}

	ImGui::End();
}

void GLErrorCheck()
{
	while (GLenum error = glGetError())
		LogMan::Log(("Caught OpenGL Error: " + std::to_string(error)), LOG_ERROR);
}

int MainUILoop()
{
	bool showDemo = true;
	while (!glfwWindowShouldClose(mainWindow))
	{
		glfwPollEvents();
		GLErrorCheck();

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
	
	glDeleteFramebuffers(0, &(viewportBuffer.fbo));
	glDeleteFramebuffers(0, &(viewportBuffer.rbo));

	TerminateMainWindow();

	return SUCCESS;
}

int StartUI()
{
	LogMan::Log("GUI Startup");
	if (!InitializeMainWindow())
	{
		//std::cout << "Failed to initialize main window" << std::endl;
		LogMan::Log("Failed to initialize main window!", LOG_ERROR);
		return FAILED_MAIN_WINDOW_INITIALIZE;
	}

	if (!InitViewport())
	{
		TerminateMainWindow();
		LogMan::Log("Failed to initialize viewport!", LOG_ERROR);
		return FAILED_VIEWPORT_CREATE;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0); //restore default framebuffer

	RecomputeWindowElementsDimensions();// mainWinWidth, mainWinHeight);
	
	LogMan::Log("GUI startup success!", LOG_SUCCESS);
	return MainUILoop();
}