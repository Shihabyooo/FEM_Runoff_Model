#pragma once
#include "GUI_Handler.hpp"
#include "GUI_Requirements.hpp"
#include "Viewport.hpp"

//const int minMainWinWidth = 800, minMainWinHeight = 600;
//const int minViewportWidth = 800, minViewportHeight = 600;

int mainWinWidth = 1280, mainWinHeight = 720;
int viewportWidth = 800, viewportHeight = 600;

GLFWwindow * mainWindow;
GLData viewportGLData;

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

Vector2 delta;
float scale;
float margin = 0.01f;
float * nodesMeshVerts = NULL;
int renderVertsCount = 3;
const double nodeDisplaySize = 0.05F;

//#pragma optimize("", off)
void UpdateNodes()
{
	std::cout << "Attemptint to update nodes with " << nodes.size() << " new nodes\n";

	if (nodes.size() < 2)
		return;

	if (nodesMeshVerts != NULL)
		delete[] nodesMeshVerts;
	renderVertsCount = nodes.size() * 9;
	nodesMeshVerts = new float[nodes.size() * 9]; //three cords * three verts per node (to draw a tri)

	delta = nodesNE - nodesSW;
	delta.x += 2.0f * margin;
	delta.y += 2.0f * margin;
	scale = Max((delta.x - 2.0f * margin) / viewPortDimensions.width, (delta.y - 2.0f * margin) / viewPortDimensions.height);

	/*std::cout << "sw: " << nodesSW.x << ", " << nodesSW.y << std::endl;
	std::cout << "ne: " << nodesNE.x << ", " << nodesNE.y << std::endl;
	std::cout << "delta: " << delta.x << ", " << delta.y << std::endl;
	std::cout << "scale: " << scale << std::endl;*/

	int counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		/*std::cout << "==================\n";
		std::cout << it->x << ", " << it->y << std::endl;*/

		Vector2 relativePos(margin + ((it->x - nodesSW.x - 0.5 * delta.x) / delta.x),
							margin + ((it->y - nodesSW.y - 0.5 * delta.x) / delta.y));
		
		//std::cout << relativePos.x << ", " << relativePos.y << std::endl;

		nodesMeshVerts[counter]		= relativePos.x;
		nodesMeshVerts[counter + 1]	= relativePos.y + nodeDisplaySize;
		nodesMeshVerts[counter + 2] = 0.0f;

		nodesMeshVerts[counter + 3] = relativePos.x + 0.5f * nodeDisplaySize;
		nodesMeshVerts[counter + 4] = relativePos.y - 0.866f * nodeDisplaySize;
		nodesMeshVerts[counter + 5] = 0.0f;

		nodesMeshVerts[counter + 6] = relativePos.x - 0.5f * nodeDisplaySize;
		nodesMeshVerts[counter + 7] = relativePos.y - 0.866f * nodeDisplaySize;
		nodesMeshVerts[counter + 8] = 0.0f;

		/*for (int i = counter; i < counter + 9; i += 3)
			std::cout << nodesMeshVerts[i]  << ", " << nodesMeshVerts[i + 1] << " , " << nodesMeshVerts[i + 2] << std::endl;*/

		counter += 9;
	}

	UpdateMesh(&viewportGLData, nodesMeshVerts, nodes.size() * 9);
}
//#pragma optimize("", on)

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
	ImGui::InputText("Mesh Nodes", meshNodes, IM_ARRAYSIZE(meshNodes));

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
		UpdateNodes();
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