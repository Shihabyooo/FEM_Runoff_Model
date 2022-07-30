#pragma once
#include "Main.hpp"
#include "LogManager.hpp"

#include "ModelInterface.hpp" //for testing only

ModelParameters ModelTestParams()
{
	ModelParameters params;

	params.demPath = "Test_Data\\DEM.tif";
	params.slopesPath = "Test_Data\\Slope_Percent.tif";
	params.fdrPath = "Test_Data\\FDR.tif";
	params.outletNode = 22;

	params.useLumpedForm = true;

	params.variablePrecipitation = false;
	params.unitTimeSeries;
	LoadTimeSeries("Test_Data\\Test_Timeseries.csv", params.unitTimeSeries);
	params.unitTimeSeries.timeUnit = TimeUnit::minute;

	params.lossModel = LossModel::scsCN;
	params.scsCN = 80;

	params.variableManningCoefficients = false;
	params.fixedManningCoeffient = 0.04; //must be positive value > 0.0
	params.manningCoefficientRasterPath = "";
	params.timeStep = 0.25; //delta T, in hours. e.g. 0.5 = 30 minutes, 1.0 = 1 hour.
	params.startTime = 0.0; //should be left at 0.0
	params.endTime = 20.0f; //hours after startTime to end simulation.

	params.useLumpedForm = true;
	params.femOmega = 0.5;
	params.solverType = Solver::Auto;
	params.residualThreshold = 0.00001;
	params.weight = 0.5;
	params.maxIterations = 5000;
	params.externalResidualTreshold = 0.000001;
	params.maxExternalIterations = 50;

	return params;
}

void GenTestMesh()
{
	MeshGeneratorParameters meshParams;
	
	meshParams.resolution = 4;
	meshParams.internalPadding = 0.01;
	meshParams.rayCastPadding = 10.0;

	GenerateMesh(meshParams);
}

int main(int argc, char ** argv)
{
	LogMan::Init(true);
	LogMan::Log("Startup.");
	
	int returnVal = 0;
	//TODO uncomment this after testing is done.
	returnVal = StartUI(1280, 720);
	
	////testing model on CLI directly
	LoadWatershedBoundary("Test_Data\\Watershed_Boundary.kml");
	GenTestMesh();
	Simulate(ModelTestParams());

	LogMan::Terminate();

	std::cin.sync(); //TODO remove
	std::cin.get(); //TODO remove

	return returnVal;
}

