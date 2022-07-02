#include "Main.hpp"
#include "LogManager.hpp"

#include "ModelInterface.hpp" //for testing only

ModelParameters ModelTestParams()
{
	ModelParameters params;

	params.demPath = "";
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
	params.fixedManningCoeffient = 0.002; //must be positive value > 0.0
	params.manningCoefficientRasterPath = "";
	params.timeStep = 0.5; //delta T, in hours. e.g. 0.5 = 30 minutes, 1.0 = 1 hour.
	params.startTime = 0.0; //should be left at 0.0
	params.endTime = 15.0f; //hours after startTime to end simulation.

	params.useLumpedForm = true;
	params.femOmega = 0.5;
	params.solverType = Solver::Auto;
	params.residualThreshold = 0.00001;
	params.weight = 1.0;
	params.maxIterations = 1000;
	params.internalResidualTreshold = 0.0001;
	params.maxInternalIterations = 100;

	return params;
}

int main(int argc, char ** argv)
{
	LogMan::Init(true);
	LogMan::Log("Startup.");
	
	int returnVal = 0;
	//TODO uncomment this after testing is done.
	//returnVal = StartUI(1280, 720);
	
	GenerateMesh("Test_Data\\Grid_Nodes.csv", 1.0);
	Simulate(ModelTestParams());

	LogMan::Terminate();

	std::cin.sync(); //TODO remove
	std::cin.get(); //TODO remove

	return returnVal;
}
