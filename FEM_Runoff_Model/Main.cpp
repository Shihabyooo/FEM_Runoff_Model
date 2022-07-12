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
	params.topographySamplingMethod = InterpolationType::linear;

	params.variablePrecipitation = false;
	params.unitTimeSeries;
	LoadTimeSeries("Test_Data\\Test_Timeseries.csv", params.unitTimeSeries);
	params.unitTimeSeries.timeUnit = TimeUnit::minute;
	params.precipitationTemporalInterpolationType = InterpolationType::linear;
	params.precipitationSpatialInterpolationType = InterpolationType::nearest;

	params.variableManningCoefficients = false;
	params.fixedManningCoeffient = 0.03; //must be positive value > 0.0
	params.manningCoefficientRasterPath = "";
	params.timeStep = 0.5; //delta T, in hours. e.g. 0.5 = 30 minutes, 1.0 = 1 hour.
	params.startTime = 0.0; //should be left at 0.0
	params.endTime = 7.5f; //hours after startTime to end simulation.

	params.useLumpedForm = true;
	params.femOmega = 0.5;
	params.solverType = Solver::Auto;
	params.residualThreshold = 0.00001;
	params.weight = 0.5;
	params.maxIterations = 1000;
	params.internalResidualTreshold = 0.0001;
	params.maxInternalIterations = 5;

	return params;
}

void TestGenerateSyntheticWatershed()
{
	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 7; j++)
			nodes.push_back(Vector2D(j * 83.33, i * 80));
	
	int counter = 0;
	for (int i = 0; i < 60; i = i + 2)
	{
		triangles.insert({ i, Triangle(i, counter, counter + 6, counter + 7, nodes) });
		triangles.insert({ i, Triangle(i, counter, counter + 1, counter + 7, nodes) });
		counter++;
	}

	for (int i = 0; i < 6; i++)
	{
		boundaryNodes.push_back(i);
		boundaryNodes.push_back(i + 36);
	}
	for (int i = 0; i < 37; i += 6)
	{
		boundaryNodes.push_back(i + 6);
		boundaryNodes.push_back(i + 11);
	}

}

int main(int argc, char ** argv)
{
	LogMan::Init(true);
	LogMan::Log("Startup.");
	
	int returnVal = 0;
	//TODO uncomment this after testing is done.
	returnVal = StartUI(1280, 720);
	

	////testing model on CLI directly
	////TestGenerateSyntheticWatershed();
	//GenerateMesh("Test_Data\\Grid_Nodes.csv", 1.0);
	//Simulate(ModelTestParams());

	LogMan::Terminate();

	//std::cin.sync(); //TODO remove
	//std::cin.get(); //TODO remove

	return returnVal;
}
