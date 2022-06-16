#pragma once
#include "GUI_Handler.hpp"
#include "GUI_Requirements.hpp"
#include "Viewport.hpp"
#include "LogManager.hpp"

int mainWinWidth = 1280, mainWinHeight = 720;

GLFWwindow * mainWindow;

ImVec4 mainBGColour = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
ImVec4 viewportBGColour = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

WindowDimensions leftPaneDimensions, logPaneDimensions, viewportDimensions;
WindowDimensions toolbarDimensions, statusBarDimensions;

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

	//disable saving of window pos (don't need it).
	io.IniFilename = NULL;
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
	//attempt to use GL 3.3
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
	glfwSetWindowSizeLimits(mainWindow, minMainWinWidth, minMainWinHeight, GLFW_DONT_CARE, GLFW_DONT_CARE); //set min screen dims
	glewInit();
	InitializeDearIMGUI(glslVersion);
	glfwSetWindowSizeCallback(mainWindow, OnMainWindowSizeChange);
	
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
	int fixedToolBarHeight = 50;
	int fixedStatusBarHeight = 30;

	statusBarDimensions = WindowDimensions(0, mainWinHeight - fixedStatusBarHeight, mainWinWidth, fixedStatusBarHeight);
	leftPaneDimensions = WindowDimensions(0, 0, fixedLeftPaneWidth, mainWinHeight - fixedStatusBarHeight);
	logPaneDimensions = WindowDimensions(leftPaneDimensions.width, mainWinHeight - fixedLogPaneHeight - fixedStatusBarHeight, mainWinWidth - fixedLeftPaneWidth, fixedLogPaneHeight);
	toolbarDimensions = WindowDimensions(fixedLeftPaneWidth, 0, mainWinWidth - fixedLeftPaneWidth, fixedToolBarHeight);
	//lastViewportSize must be set before we update viewportDimensions.
	lastViewportSize.x = viewportDimensions.width;
	lastViewportSize.y = viewportDimensions.height;
	viewportDimensions = WindowDimensions(leftPaneDimensions.width, fixedToolBarHeight, mainWinWidth - fixedLeftPaneWidth, mainWinHeight - fixedLogPaneHeight - fixedStatusBarHeight - fixedToolBarHeight);

	//update viewport renderer to use new dimensions.
	UpdateViewport();
}

//TODO move window/viewport rendering/painting function to their own header/source files.
//char meshNodes[260] = "Path to grid nodes file.";
char meshNodes[260] = "Test_Nodes_R1.csv";
char demFilePath[260] = "Test_Raster.tif";

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
	ImGui::PushItemWidth(-1);
	ImGui::InputText("DEM File Path", demFilePath, IM_ARRAYSIZE(demFilePath));
	ImGui::PopItemWidth();

	if (ImGui::Button("Browse for DEM directory"))
	{

	}

	if (ImGui::Button("Load DEM.", ImVec2(100, 50)))
	{
		TestLoadDEM(demFilePath);
	}
		
	//Other input
	ImGui::NewLine();
	if (ImGui::Button("Run Simulation!", ImVec2(100, 50)))
	{
		TestSimulate(meshNodes);
		SetViewBounds(nodesSW, nodesNE);
		UpdateViewport();
	}

	ImGui::End();
}

void DrawToolbar()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	ImGui::SetNextWindowPos(ImVec2(toolbarDimensions.positionX, toolbarDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(toolbarDimensions.width, toolbarDimensions.height), ImGuiCond_Always);

	ImGui::Begin("Toolbar", NULL, windowFlags);
	

	ImGui::End();
}

void DrawStatusBar()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	ImGui::SetNextWindowPos(ImVec2(statusBarDimensions.positionX, statusBarDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(statusBarDimensions.width, statusBarDimensions.height), ImGuiCond_Always);

	ImGui::Begin("Statusbar", NULL, windowFlags);
	
	if (isHoveringViewport)
		ImGui::Text("Hover Pos: |%g, %g|", currenViewportHoverPos.x, currenViewportHoverPos.y);
	else
		ImGui::Text("Hover Pos: |--, --|");
	ImGui::SameLine(300);
	ImGui::Text("View Rect: |%g, %g : %g, %g|", viewBounds[0].x, viewBounds[0].y, viewBounds[1].x, viewBounds[1].y);
	ImGui::SameLine();
	ImGui::Text("Scale: %.4g", scale);
	ImGui::SameLine();
	ImGui::Text("Screen Aspect: %.4g", screenAspectRatio);
	ImGui::SameLine();
	ImGui::Text("World Aspect: %.4g", worldAspectRatio);
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
		
		//update viewport to current main window size
		glViewport(0, 0, mainWinWidth, mainWinHeight); 

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		glClearColor(mainBGColour.x * mainBGColour.w, mainBGColour.y * mainBGColour.w, mainBGColour.z * mainBGColour.w, mainBGColour.w);
		glClear(GL_COLOR_BUFFER_BIT);

		DrawLeftPane();
		DrawLogPane();
		DrawToolbar();
		DrawStatusBar();
		DrawViewport();

		ImGui::ShowDemoWindow(&showDemo);

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
	return StartUI(minMainWinWidth, minMainWinHeight);
}

int StartUI(unsigned int const startResX, unsigned int const startResY)
{
	LogMan::Log("GUI Startup");
	mainWinWidth = startResX;
	mainWinHeight = startResY;

	if (!InitializeMainWindow())
	{
		LogMan::Log("Failed to initialize main window!", LOG_ERROR);
		return FAILED_MAIN_WINDOW_INITIALIZE;
	}

	if (!InitViewport())
	{
		TerminateMainWindow();
		LogMan::Log("Failed to initialize viewport!", LOG_ERROR);
		return FAILED_VIEWPORT_CREATE;
	}

	//update main window dimensions then compute GUI elements positions/dimensions.
	//We have to do this here because RecomputeWindowElementsDimensions() calls UpdateViewport(), which assumes viewport is initialized.
	glfwGetWindowSize(mainWindow, &mainWinWidth, &mainWinHeight);
	RecomputeWindowElementsDimensions();

	glBindFramebuffer(GL_FRAMEBUFFER, 0); //restore default framebuffer
	
	LogMan::Log("GUI startup success!", LOG_SUCCESS);
	return MainUILoop();
}

