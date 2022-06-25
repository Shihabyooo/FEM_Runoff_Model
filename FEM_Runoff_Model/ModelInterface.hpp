#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "Globals.hpp"
#include "DelauneyTriangulation.hpp"



//Process is as follows:
//Recieve a model parameters from GUI (or CLI)
//Load required rasters and datasets that need loading.
//Data check pass
//Actual model run
//Compute fixed global matrices
	//Capacitance Matrix (C)
	//X-Conductance Matrix (Psi-X)
	//Y-Conductance Matrix (Psi-Y)
	//Beta vector
//Set initial heads to zero (Dry conditions).
//Loop from start time to end time
	//Compute effective rainfall for each element
	//Compute qx and qy vectors using previous pass heads (for first loop, it's the initial head)
	//Internal loop
		//Solve the system of equations.
		//Check the resulting heads with the ones assumed for qx and qy, if the difference is large, recompute qx and qy, loop again.
		//if difference is accepable, break internal loop.
	//handle data storage for results, residuals and any relative statistics, prepare for next loop.
//Display results and return control to user.

extern std::unordered_map<int, Triangle> triangles;
extern std::vector<Vector2> nodes;
extern std::vector<int> boundaryNodes;
extern Vector2 nodesSW, nodesNE;

bool GenerateMesh(std::string const & nodesPath);
bool Simulate(ModelParameters const & params);
