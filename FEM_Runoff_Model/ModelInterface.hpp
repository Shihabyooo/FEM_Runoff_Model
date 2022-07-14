#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "Globals.hpp"
#include "MeshGenerator.hpp"


//TODO the current time approximation uses a double for the omega variable. While it is mathematically correct, practically, the values\
should typically either be 0.0, 0.5 or 1.0. In the first and last case, the RHS computation can be significantly sped up by omitting\
that would -in each case- result in a vector of zeroes (see approximate solution formulation). It may make sense hence to instead use a\
"mode" variable (i.e. an int = 0, 1 or 2) instead of a double, and a switch statement to do an optimized RHS computation based on the mode\
e.g. if mode = 0 (backwards difference), qx, qy and Pe at t+dt are unnecessary, for mode = 2 (forward diff), qx, qy and Pe at t are uncessary.

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

//Note: Rainfall is expected in mm, but the spatial units are meters (squared), so divide Pe by 1000 before multiplying with beta.



//extern ElementsList elements;
extern std::unordered_map<size_t, Rectangle> rectangles;
extern std::unordered_map<size_t, Triangle> triangles;
extern Triangle superTriangles[2];
extern std::vector<Vector2D> nodes;
extern std::vector<size_t> boundaryNodes;
extern Vector2D nodesSW, nodesNE;
extern Vector2D shedSW, shedNE;
//extern size_t exitNode;
extern std::vector<Vector2D> shedBoundary;


bool LoadWatershedBoundary(std::string const & boundaryPath);
bool GenerateMesh(std::string const & nodesPath, double superTrianglePadding);
bool GenerateGridMesh(size_t resolution, double internalPadding, double raycastPadding);
Triangle const * GetElementContainingPoint(Vector2D const & pos);
bool UpdateNode(size_t id, Vector2D const & newPos);

bool LoadTimeSeries(std::string const & path, TimeSeries & ts);
bool Simulate(ModelParameters const & params);
