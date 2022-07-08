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

int minLeftPaneWidth = 250;
int maxLeftPaneWidth = 1000;
int minLogPaneHeight = 200;
int maxLogPaneHeight = 500;

//int fixedLeftPaneWidth = 250;
//int fixedLogPaneHeight = 200;
int fixedToolBarHeight = 60;
int fixedStatusBarHeight = 35;

Icon iconNodeMove, iconNodeView, iconElemView;

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

	io.IniFilename = "gui_conf.ini";
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

bool LoadIcons()
{
	iconNodeMove = std::move(Icon(FileIO::LoadImage("Icons\\Icon_moveNode.png")));
	iconNodeView = std::move(Icon(FileIO::LoadImage("Icons\\Icon_viewNode.png")));
	iconElemView = std::move(Icon(FileIO::LoadImage("Icons\\Icon_viewElement.png")));

	return iconNodeMove.isInit && iconNodeView.isInit && iconElemView.isInit;
}

void TerminateMainWindow()
{
	glfwDestroyWindow(mainWindow);
	glfwTerminate();
}

void RecomputeWindowElementsDimensions()// int newMainWinWidth, int newMainWinHeight)
{
	//When minimized, width and height are zero. This can cause issue with viewport renderer. So we simply skip updating dimensions
	if (mainWinHeight == 0 || mainWinWidth == 0) 
		return;	

	//Set fixed width for left pane.
	//Set fixed height for log pane.
	//left pane height = mainWindow height.
	//log pane width = mainwindow width - fixed left pane width.
	//viewport takes remainin area.
	
	statusBarDimensions = WindowDimensions(0, mainWinHeight - fixedStatusBarHeight, mainWinWidth, fixedStatusBarHeight);
	leftPaneDimensions = WindowDimensions(0, 0, leftPaneDimensions.width, mainWinHeight - fixedStatusBarHeight);
	logPaneDimensions = WindowDimensions(leftPaneDimensions.width, mainWinHeight - logPaneDimensions.height - fixedStatusBarHeight, mainWinWidth - leftPaneDimensions.width, logPaneDimensions.height);
	toolbarDimensions = WindowDimensions(leftPaneDimensions.width, 0, mainWinWidth - leftPaneDimensions.width, fixedToolBarHeight);
	
	lastViewportSize.x = viewportDimensions.width;
	lastViewportSize.y = viewportDimensions.height;
	viewportDimensions = WindowDimensions(leftPaneDimensions.width,
		fixedToolBarHeight,
		mainWinWidth - leftPaneDimensions.width,
		mainWinHeight - logPaneDimensions.height - fixedStatusBarHeight - fixedToolBarHeight);

	//update viewport renderer to use new dimensions.
	UpdateViewport();
}


//char meshNodes[260] = "Mesh Nodes Coordinates File Path";
//char meshNodes[260] = "Test_Nodes_G.csv";
char meshNodes[260] = "Test_Data\\Grid_Nodes.csv";
char demFilePath[260] = "DEM Raster Path";
//char slopeFilePath[260] = "Slopes Raster Path";
char slopeFilePath[260] = "Test_Data\\Slope_Percent.tif";
//char fdrFilePath[260] = "Flow Direction Raster Path";
char fdrFilePath[260] = "Test_Data\\FDR.tif";

char precipRasterDir[260] = "Precipitation Rasters Directory Path";
//char percipTSPath[260] = "Precipitation Time-Series File Path";
char percipTSPath[260] = "Test_Data\\Test_Timeseries.csv";
char manningFilePath[260] = "Manning Raster File Path";

char superTriPadding[24] = "1.0";

int timeSeriesSize = 3;
TimeSeries inputTS;

bool useFixedManning = true;
char fixedManningCoef[12] = "0.002";

char startTime[24] = "0.0";
char endTime[24] = "15";
char deltaTime[24] = "0.5";

char femOmega[12] = "0.5";
bool useLumped = true;

const char* solvers[] = { "Auto", "Simple", "Gaussian", "Jacobi", "SOR", "PCG", "BiCG", "CGS" };// , "GMRES" };
const char* precipitationInput[] = {"Single Time-Series", "Gridded Time-Series" };
const char* interpMethods1D[] = {"Nearest", "Linear", "Cubic"};
const char* interpMethods2D[] = { "Nearest", "Bilinear", "Bicubic" };
const char* timeUnits[] = { "Seconds", "Minutes", "Hours", "Days" };
int selectedPrecipInput = 0;
int selectedPrecipTempoInterp = 1;
int selectedPrecipSpaceInterp = 1;
int selectedTopoInterp = 1;
int selectedTimeUnit = 1;
int selectedSolver = 0;

char solverResidual[12] = "0.0001";
char solverWeight[12] = "0.5";
int solverMaxIteration = 1000;
char internalResidual[12] = "0.0001";
int internalMaxIteration = 100;


void FillParametersStruct(ModelParameters & params)
{
	LogMan::Log("Processing Input Parameters");

	params.demPath = demFilePath;
	params.slopesPath = slopeFilePath;
	params.fdrPath = fdrFilePath;
	params.topographySamplingMethod = static_cast<InterpolationType>(selectedTopoInterp);
	params.variablePrecipitation = selectedPrecipInput == 1;

	//params.unitTimeSeries = ; 
	params.precipitationTemporalInterpolationType = static_cast<InterpolationType>(selectedPrecipTempoInterp);
	params.precipitationSpatialInterpolationType = static_cast<InterpolationType>(selectedPrecipSpaceInterp);
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
	params.internalResidualTreshold = atof(internalResidual);
	params.maxInternalIterations = internalMaxIteration;

	params.unitTimeSeries = inputTS;
	params.unitTimeSeries.timeUnit = static_cast<TimeUnit>(selectedTimeUnit);
	
	LogMan::Log("Finished processing Input Parameters", LOG_SUCCESS);
}

void DrawLeftPane()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoMove;
	
	ImGui::SetNextWindowSizeConstraints(ImVec2(minLeftPaneWidth, -1), ImVec2(maxLeftPaneWidth, -1));

	ImGui::Begin("Input", NULL, windowFlags);

	if (leftPaneDimensions.UpdateWidth(ImGui::GetWindowSize().x))
		RecomputeWindowElementsDimensions();

	ImGui::SetWindowPos(ImVec2(leftPaneDimensions.positionX, leftPaneDimensions.positionY), ImGuiCond_Always);
	ImGui::SetWindowSize(ImVec2(leftPaneDimensions.width, leftPaneDimensions.height), ImGuiCond_Always);

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

		if (ImGui::Button("Browse##nodes"))
		{
			//TODO spawn file browser here
		}
		
		ImGui::Text("SuperTriangle padding");
		ImGui::SameLine();
		ImGui::InputText("##superTriPad", superTriPadding, 24, ImGuiInputTextFlags_CharsDecimal);

		if (ImGui::Button("Generate Mesh", ImVec2(100, 50)))
		{
			GenerateMesh(meshNodes, atof(superTriPadding));
			SetViewBounds(nodesSW, nodesNE);
			UpdateViewport();
		}
		
		ImGui::Text("DEM");
		ImGui::InputText("##demFilePath", demFilePath, IM_ARRAYSIZE(demFilePath));
		
		if (ImGui::Button("Browse##dem"))
			LogMan::Log("Not yet impltemented!", LOG_WARN);

		ImGui::Text("Slopes");
		ImGui::InputText("##slopesFilePath", slopeFilePath, IM_ARRAYSIZE(slopeFilePath));

		if (ImGui::Button("Browse##slopes"))
			LogMan::Log("Not yet impltemented!", LOG_WARN);

		ImGui::Text("Flow Direction");
		ImGui::InputText("##fdrFilePath", fdrFilePath, IM_ARRAYSIZE(fdrFilePath));

		if (ImGui::Button("Browse##fdr"))
			LogMan::Log("Not yet impltemented!", LOG_WARN);

		ImGui::Text("Spatial interpolation method");
		ImGui::Combo("##InterpT2D", &selectedTopoInterp, interpMethods2D, IM_ARRAYSIZE(interpMethods2D));
		ImGui::PopItemWidth();
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Precipitation"))
	{
		ImGui::PushItemWidth(-1);
		ImGui::Text("Precipitation input method");		
		ImGui::Combo("##precipInput", &selectedPrecipInput, precipitationInput, IM_ARRAYSIZE(precipitationInput));

		if (selectedPrecipInput == 0)
		{
			ImGui::InputInt("Series length", &timeSeriesSize);
			timeSeriesSize = Max(timeSeriesSize, 2);
			ImGui::Text("Time-series temporal interpolation method");
			ImGui::Combo("##InterpP1D", &selectedPrecipTempoInterp, interpMethods1D, IM_ARRAYSIZE(interpMethods1D));

			ImGui::Text("Time-series time units");
			ImGui::Combo("##tsTimeUnits", &selectedTimeUnit, timeUnits, IM_ARRAYSIZE(timeUnits));

			ImGui::Text("Time-series file path");
			ImGui::InputText("##tsFilePath", percipTSPath, IM_ARRAYSIZE(percipTSPath));			
	
			if (ImGui::Button("Load Time-Series"))
			{
				LoadTimeSeries(percipTSPath, inputTS);
				timeSeriesSize = inputTS.size;
			}

			if (ImGui::CollapsingHeader("Time-Series", ImGuiTreeNodeFlags_None))
			{
				inputTS.AdjustSize(timeSeriesSize); //will do nothing if size is same
			
				ImGui::BeginTable("##timeSeriesTable", 3, ImGuiTableFlags_Resizable, ImVec2(leftPaneDimensions.width - 10, 100));
				
				ImGui::TableNextColumn();
				ImGui::Text(" ");
				ImGui::TableNextColumn();
				std::string timeColumnHeader = "Time (" + std::string(timeUnits[selectedTimeUnit]) + ")";
				ImGui::Text(timeColumnHeader.c_str());
				ImGui::TableNextColumn();
				ImGui::Text("Precipitation (mm)");

				ImGui::TableNextColumn();
				ImGui::Text("1");
				ImGui::TableNextColumn();
				ImGui::Text("0");
				ImGui::TableNextColumn();
				ImGui::Text("0.0");

				for (int i = 1; i < timeSeriesSize; i++)
				{
					ImGui::TableNextColumn();
					ImGui::Text(std::to_string(i + 1).c_str());

					ImGui::TableNextColumn();

					std::pair<int, double> tempPair(inputTS.series[i].first, inputTS.series[i].second); //to force input to positive values.

					ImGui::InputInt(std::string("##precipTime" + std::to_string(i)).c_str(), &tempPair.first, 0);
					ImGui::TableNextColumn();
					ImGui::InputDouble(std::string("##precipVal" + std::to_string(i)).c_str(), &tempPair.second);
					
					//force input to positive values
					inputTS.series[i].first = Max(tempPair.first, 0);
					inputTS.series[i].second = Max(tempPair.second, 0.0);
				}
				ImGui::EndTable();
			}
		}
		else if (selectedPrecipInput == 1)
		{
			ImGui::Text("Precipitation rasters directories");
			ImGui::InputText("##precipRasDir", precipRasterDir, IM_ARRAYSIZE(precipRasterDir));
			ImGui::Text("Time-series temporal interpolation method");
			ImGui::Combo("##InterpP1D", &selectedPrecipTempoInterp, interpMethods1D, IM_ARRAYSIZE(interpMethods1D));

			ImGui::Text("Precipitation spatial interpolation method");
			ImGui::Combo("##InterpP2D", &selectedPrecipSpaceInterp, interpMethods2D, IM_ARRAYSIZE(interpMethods2D));
		}
		ImGui::PopItemWidth();
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Hydraulics"))
	{
		ImGui::PushItemWidth(-1);
		ImGui::Checkbox("Use Fixed Manning Coefficient", &useFixedManning);
		if (useFixedManning)
		{
			ImGui::Text("Manning Coefficient:");
			ImGui::SameLine();
			ImGui::InputText("##fixedManningCoef", fixedManningCoef, 12, ImGuiInputTextFlags_CharsDecimal);
		}
		else
		{
			ImGui::Text("Manning Coefficients Raster");
			ImGui::InputText("##ManningRaster", manningFilePath, IM_ARRAYSIZE(manningFilePath));
		}
		ImGui::PopItemWidth();
	}

	ImGui::NewLine();
	ImGui::Separator();

	if (ImGui::CollapsingHeader("Simulation Run Parameters"))
	{
		ImGui::PushItemWidth(-1);
		ImGui::Text("Start time:");
		ImGui::SameLine();
		ImGui::InputText("##startTime", startTime, 24, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("End time:");
		ImGui::SameLine();
		ImGui::InputText("##endTime", endTime, 24, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("Time step:");
		ImGui::SameLine();
		ImGui::InputText("##deltaTime", deltaTime, 24, ImGuiInputTextFlags_CharsDecimal);
		ImGui::NewLine();
		
		ImGui::Checkbox("Use Lumped formulation", &useLumped);
		ImGui::Text("Omega:");
		ImGui::SameLine();
		ImGui::InputText("##femOmega", femOmega, 12, ImGuiInputTextFlags_CharsDecimal);
		
		ImGui::NewLine();
		ImGui::Text("Solver:");
		ImGui::Combo("##Solver", &selectedSolver, solvers, IM_ARRAYSIZE(solvers));

		ImGui::Text("Max residual:");
		ImGui::SameLine();
		ImGui::InputText("##residual", solverResidual, 12, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("Max iterations:");
		ImGui::SameLine();
		ImGui::InputInt("##iterations", &solverMaxIteration);
	
		ImGui::Text("Max residual (Internal):");
		ImGui::SameLine();
		ImGui::InputText("##residualInternal", internalResidual, 12, ImGuiInputTextFlags_CharsDecimal);
		ImGui::Text("Max iterations (Internal):");
		ImGui::SameLine();
		ImGui::InputInt("##iterationsInteral", &internalMaxIteration);



		if (selectedSolver == 0 || selectedSolver == 3 || selectedSolver == 4)
		{
			ImGui::Text("Weight:");
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

ImVec2 buttonSize(32, 32);
ImVec4 buttonBGColour(1.0, 1.0f, 1.0f, 0.0f);
ImVec4 buttonTintColour(1.0f, 1.0f, 1.0f, 1.0f);
ImVec4 activeButtonTintColour(0.45f, 0.75f, 1.0f, 1.0f);

void DrawButton(Icon const & icon, ToolMode toolMode)
{
	if (ImGui::ImageButton((void*)icon.texPtr,
							buttonSize,
							ImVec2(0, 0),
							ImVec2(1, 1),
							-1,
							buttonBGColour,
							GetActiveTool() == toolMode ? activeButtonTintColour : buttonTintColour))
	{
		SwitchActiveTool(toolMode);
	}
}

void DrawToolbar()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	
	ImGui::Begin("Toolbar", NULL, windowFlags);
	
	ImGui::SetWindowPos(ImVec2(toolbarDimensions.positionX, toolbarDimensions.positionY), ImGuiCond_Always);
	ImGui::SetWindowSize(ImVec2(toolbarDimensions.width, toolbarDimensions.height), ImGuiCond_Always);
	
	//move nodes
	DrawButton(iconNodeMove, ToolMode::MoveNode);
	
	//view nodes
	ImGui::SameLine();
	DrawButton(iconNodeView, ToolMode::ViewNode);
	
	//view elements
	ImGui::SameLine();
	DrawButton(iconElemView, ToolMode::ViewElement);
	
	//ImGui::SameLine();
	//if (ImGui::ImageButton(NULL, buttonSize, ImVec2(0, 0), ImVec2(1, 1), -1, buttonBGColour, buttonTintColour))
	//{

	//}

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
	
	if (!LoadIcons())
		LogMan::Log("Warning! There was an issue loading icons.", LOG_WARN);
	
	LogMan::Log("GUI startup success!", LOG_SUCCESS);
	return MainUILoop();
}

