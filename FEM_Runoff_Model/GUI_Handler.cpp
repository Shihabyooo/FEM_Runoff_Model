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

int fixedLeftPaneWidth = 250;
int fixedLogPaneHeight = 200;
int fixedToolBarHeight = 50;
int fixedStatusBarHeight = 30;

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


char meshNodes[260] = "Mesh Nodes Coordinates File Path";
char demFilePath[260] = "DEM Raster Path";
char precipRasterDir[260] = "Precipitation Rasters Directory Path";
char manningFilePath[260] = "Manning Raster File Path";

char fixedPrecipVal[24] = "";
int timeSeriesSize = 3;

bool useFixedManning = true;
char fixedManningCoef[12] = "";

char startTime[24] = "0.0";
char endTime[24] = "";
char deltaTime[24] = "";

char femOmega[12] = "0.5";
bool useLumped = true;

const char* solvers[] = { "Auto", "Simple", "Gaussian", "Jacobi", "SOR", "PCG", "BiCG", "CGS" };// , "GMRES" };
const char* precipitationInput[] = { "Single, fixed value", "Single Time-Series", "Gridded Time-Series" };
const char* interpMethods1D[] = {"Nearest", "Linear", "Cubic"};
const char* interpMethods2D[] = { "Nearest", "Bilinear", "Bicubic" };
int selectedPrecipInput = 0;
int selectedPrecipTempoInterp = 1;
int selectedPrecipSpaceInterp = 1;
int selectedTopoInterp = 1;
int selectedSolver = 0;

char solverResidual[12] = "0.01";
char solverWeight[12] = "1.0";
int solverMaxIteration = 1000;


void FillParametersStruct(ModelParameters & params)
{
	params.demPath = demFilePath;
	params.topographySamplingMethod = static_cast<InterpolationType>(selectedTopoInterp);
	params.variablePrecipitation = selectedPrecipInput == 0;

	//params.unitTimeSeries = ; 
	params.precipitationTemporalInterpolationType = static_cast<InterpolationType>(selectedPrecipTempoInterp);
	params.precipitationSpatialInterpolationType = static_cast<InterpolationType>(selectedPrecipSpaceInterp);
	params.fixedPrecipitationValue = atof(fixedPrecipVal);
	params.variableManningCoefficients = !useFixedManning;
	params.fixedManningCoeffient = atof(fixedManningCoef);
	params.manningCoefficientRasterPath = manningFilePath;

	//params.useBuiltInLossModel = ;
	//params.useHydrologicClassGrid = ; //if true, hydrologic class raster must be set
	//params.hydrologicClassRaster = ;
	params.timeStep = atof(deltaTime);
	params.startTime = atof(startTime);
	params.endTime = atof(endTime);

	params.useLumpedForm = useLumped;
	params.femOmega = atof(femOmega);

	params.solverType = static_cast<Solver>(selectedSolver);
	params.residualThreshold = atof(solverResidual);
	params.weight = atof(solverWeight);
	params.maxIterations = solverMaxIteration < 0 ? 0 : solverMaxIteration;

}

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

	if (ImGui::CollapsingHeader("Topography and Geometry"))
	{
		//Geometry data
		ImGui::Text("Mesh Nodes");
		ImGui::PushItemWidth(-1);
		ImGui::InputText("Mesh Nodes", meshNodes, IM_ARRAYSIZE(meshNodes));
		ImGui::PopItemWidth();
		if (ImGui::Button("Browse for geometry directory"))
		{
			//TODO spawn file browser here
		}
		
		ImGui::NewLine();

		if (ImGui::Button("Generate Mesh", ImVec2(100, 50)))
		{
			GenerateMesh(meshNodes);
			SetViewBounds(nodesSW, nodesNE);
			UpdateViewport();
		}

		/*DrawFileBrowser();
		DrawFileList(geometryFilePath, &geometryNames, &selectedGeometry, DataType::geometry, false, GEOMETRY_LIST_ID);
		ImGui::NewLine();*/

		ImGui::Text("DEM");
		ImGui::PushItemWidth(-1);
		ImGui::InputText("##demFilePath", demFilePath, IM_ARRAYSIZE(demFilePath));
		ImGui::PopItemWidth();

		if (ImGui::Button("Browse for DEM directory"))
			LogMan::Log("Not yet impltemented!", LOG_WARN);

		ImGui::Text("Spatial interpolation method");
		ImGui::Combo("##InterpT2D", &selectedTopoInterp, interpMethods2D, IM_ARRAYSIZE(interpMethods2D));
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Precipitation"))
	{
		ImGui::Text("Precipitation input method"); //TODO text warpping
		ImGui::Combo("##precipInput", &selectedPrecipInput, precipitationInput, IM_ARRAYSIZE(precipitationInput));

		if (selectedPrecipInput == 0)
		{
			ImGui::Text("Fixed Inrecemental Precipitation");
			ImGui::InputText("##fixedPercipVal", fixedPrecipVal, 24, ImGuiInputTextFlags_CharsDecimal);
		}
		else if (selectedPrecipInput == 1)
		{
			ImGui::InputInt("Series length", &timeSeriesSize);
			ImGui::Text("Time-series temporal interpolation method"); //TODO text warpping
			ImGui::Combo("##InterpP1D", &selectedPrecipTempoInterp, interpMethods1D, IM_ARRAYSIZE(interpMethods1D));
			if (ImGui::CollapsingHeader("Time-Series", ImGuiTreeNodeFlags_None))
			{
				//TODO time series input (HEC-HMS-esque)
			}
		}
		else if (selectedPrecipInput == 2)
		{
			ImGui::Text("Precipitation rasters directories");
			ImGui::PushItemWidth(-1);
			ImGui::InputText("##precipRasDir", precipRasterDir, IM_ARRAYSIZE(precipRasterDir));
			ImGui::PopItemWidth();
			ImGui::Text("Time-series temporal interpolation method"); //TODO text warpping
			ImGui::Combo("##InterpP1D", &selectedPrecipTempoInterp, interpMethods1D, IM_ARRAYSIZE(interpMethods1D));

			ImGui::Text("Precipitation spatial interpolation method"); //TODO text warpping
			ImGui::Combo("##InterpP2D", &selectedPrecipSpaceInterp, interpMethods2D, IM_ARRAYSIZE(interpMethods2D));
		}
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Hydraulics"))
	{
		ImGui::Checkbox("Use Fixed Manning Coefficient", &useFixedManning);
		if (useFixedManning)
		{
			ImGui::PushItemWidth(-1);
			ImGui::Text("Manning Coefficient");
			ImGui::SameLine();
			ImGui::InputText("##fixedManningCoef", fixedManningCoef, 12, ImGuiInputTextFlags_CharsDecimal);
			ImGui::PopItemWidth();
		}
		else
		{
			ImGui::Text("Manning Coefficients Raster");
			ImGui::PushItemWidth(-1);
			ImGui::InputText("##ManningRaster", manningFilePath, IM_ARRAYSIZE(manningFilePath));
			ImGui::PopItemWidth();
		}
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Simulation Run Parameters"))
	{
		ImGui::PushItemWidth(-1);
		ImGui::Text("Start time");
		ImGui::SameLine();
		ImGui::InputText("##startTime", startTime, 24, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("End time");
		ImGui::SameLine();
		ImGui::InputText("##endTime", endTime, 24, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("Time step");
		ImGui::SameLine();
		ImGui::InputText("##deltaTime", deltaTime, 24, ImGuiInputTextFlags_CharsDecimal);
		ImGui::NewLine();
		
		ImGui::Checkbox("Use Lumped formulation", &useLumped);
		ImGui::Text("Time difference scheme (Omega value)");
		ImGui::SameLine();
		ImGui::InputText("##femOmega", femOmega, 12, ImGuiInputTextFlags_CharsDecimal);

		ImGui::PopItemWidth();
		
		ImGui::NewLine();
		ImGui::Text("Solver");
		ImGui::Combo("##Solver", &selectedSolver, solvers, IM_ARRAYSIZE(solvers));

		ImGui::PushItemWidth(-1);
		ImGui::Text("Max residual");
		ImGui::SameLine();
		ImGui::InputText("##residual", solverResidual, 12, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("Iterations");
		ImGui::SameLine();
		ImGui::InputInt("##iterations", &solverMaxIteration);
		

		if (selectedSolver == 3 || selectedSolver == 4)
		{
			ImGui::Text("Weight");
			ImGui::SameLine();
			ImGui::InputText("##weight", solverWeight, 12, ImGuiInputTextFlags_CharsDecimal);
		}
		ImGui::PopItemWidth();
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::Button("Run Simulation!", ImVec2(100, 50)))
	{
		//create params from current input then pass to Simulate()
		ModelParameters newParams;
		FillParametersStruct(newParams);
		Simulate(newParams);
	}

	ImGui::End();
}

void DrawToolbar()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	ImGui::SetNextWindowPos(ImVec2(toolbarDimensions.positionX, toolbarDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(toolbarDimensions.width, toolbarDimensions.height), ImGuiCond_Always);

	ImGui::Begin("Toolbar", NULL, windowFlags);
	ImVec2 buttonSize(32, 32);
	ImVec4 buttongBGColour(1.0, 1.0f, 1.0f, 1.0f);
	ImVec4 buttonTintColour(1.0f, 0.0f, 0.0f, 1.0f);
	
	ImGuiIO& io = ImGui::GetIO(); //test
	ImTextureID my_tex_id = io.Fonts->TexID; //test

	if (ImGui::ImageButton(my_tex_id, buttonSize, ImVec2(0, 0), ImVec2(1, 1), -1, buttongBGColour, buttonTintColour))
	{

	}
	ImGui::SameLine();
	if (ImGui::ImageButton(my_tex_id, buttonSize, ImVec2(0, 0), ImVec2(1, 1), -1, buttongBGColour, buttonTintColour))
	{

	}
	ImGui::SameLine();
	if (ImGui::ImageButton(my_tex_id, buttonSize, ImVec2(0, 0), ImVec2(1, 1), -1, buttongBGColour, buttonTintColour))
	{

	}
	ImGui::SameLine();
	if (ImGui::ImageButton(my_tex_id, buttonSize, ImVec2(0, 0), ImVec2(1, 1), -1, buttongBGColour, buttonTintColour))
	{

	}

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

