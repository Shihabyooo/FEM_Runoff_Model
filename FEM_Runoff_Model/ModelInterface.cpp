#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//std::unordered_map<int, Triangle> triangles;
//std::vector<Vector2> nodes;
std::unordered_map<int, Triangle> triangles;
std::vector<Vector2> nodes;
std::vector<int> boundaryNodes;
Vector2 nodesSW, nodesNE;

//Matrices and Vectors
Matrix_f32 globalC, globalPsiX, globalPsiY;
int demID, manningRasterID;

void TestBoundingBox()
{
	//std::cout << "\n computing bounding box in model interface\n";
	nodesSW = nodes[0];
	nodesNE = nodes[0];

	for (auto it = nodes.begin(); it < nodes.end(); it++)
	{
		Print(*it);
		nodesSW.x = Min(it->x, nodesSW.x);
		nodesSW.y = Min(it->y, nodesSW.y);
		nodesNE.x = Max(it->x, nodesNE.x);
		nodesNE.y = Max(it->y, nodesNE.y);
	}

	LogMan::Log("Loaded nodes with bounds: "
				+ std::to_string(nodesSW.x) + ", " + std::to_string(nodesSW.y) + " and "
				+ std::to_string(nodesNE.x) + ", " + std::to_string(nodesNE.y));
}

void ConstructGlobalConductanceMatrices(double deltaT)
{
	//Psi-X and Psi-Y (for each element) are 3x3 matrices.
	//[Psi-X_e] = 1/6 * delta T *	|	yj-yk	yk-yi	yi-yk	|
	//								|	yj-yk	yk-yi	yi-yk	|
	//								|	yj-yk	yk-yi	yi-yk	|

	//[Psi-Y_e] = 1/6 * delta T *	|	xk-xj	xi-xk	xk-xi	|
	//								|	xk-xj	xi-xk	xk-xi	|
	//								|	xk-xj	xi-xk	xk-xi	|
	//Where x(i/j/k) and y(i/j/k) are the x, y coord of the i/j/kth node.

	//Create global  Psi-X and Psi-y of size (node x nodes), zero initial value.
	//loop over each triangle in mesh
		//loop over permutations of vertIDs (3 verts per triangle = 9 permutations)
			//Psi-X[permutation] += Psi-X_e[localized permutation], e.g. if triangle is verts 3, 5, 8  (i, j, k)\
			Psi-x[3][3] += Psi-x_e[0][0] = 1/6 * dT * (yj-yk), and\
			Psi-x[5][8] += Psi-x_e[1][2] = 1/6 * dT * (yK-yI), etc
			//Ditto for Psi-Y

}

void ConstructGlobalCapacitanceMatrix()
{
	//Capacitance matrix for each element is 3x3 matrix
	//[C_e] = A/3 * |	1	0	0	|
	//				|	0	1	0	|
	//				|	0	0	1	|
	//Where A is the area of element.

	//Construct global matrix similar to procedure in ConstructGlobalConductanceMatrices()
}

void ConstructGlobalBetaMatrix()
{
	//TODO figure out how to construct this matrix (or vector?)

	//Beta matrix for each element is a 3x1 vector
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.

	////Construct global matrix similar to proceudre in ConstructGlobalConductanceMatrices()
	//Nope. This isn't correct...	
}

bool CheckParameters(ModelParameters const & params)
{
	bool status = true;

	//Check fails

	//TODO for raster file checks, check also that they are georeffed rasters and inclusive of the mesh boundary.
	if (triangles.size() < 1) //existence of tris implies existence of nodes (though I probably need to protect both from change outside this file)
	{
		LogMan::Log("ERROR! No loaded mesh.", LOG_ERROR);
		status = false;
	}

	if (!FileExists(params.demPath))
	{
		LogMan::Log("ERROR! Must supply a terrain DEM for the region.", LOG_ERROR);
		status = false;
	}

	if (params.variablePrecipitation && params.unitTimeSeries == NULL) //TODO check for time series rasters also goes here
	{
		LogMan::Log("ERROR! Must set a time series for preciptation when using variable precipitation", LOG_ERROR);
		status = false;
	}
	else if (!params.variablePrecipitation && params.fixedPrecipitationValue <= 0.0)
	{
		LogMan::Log("ERROR! Must set a positive precipitation value when using fixed precipitation", LOG_ERROR);
		status = false;
	}

	if (params.variableManningCoefficients && !FileExists(params.manningCoefficientRasterPath))
	{
		LogMan::Log("ERROR! Must set a Manning coefficients raster when using variable manning coefficients", LOG_ERROR);
		status = false;
	}
	else if (!params.variableManningCoefficients && params.fixedManningCoeffient <= 0.0)
	{
		LogMan::Log("ERROR! Must set a Manning coefficient  when using fixed manning coefficients", LOG_ERROR);
		status = false;
	}

	if (params.timeStep <= 0.0)
	{
		LogMan::Log("ERROR! Time step must be a positive, real number", LOG_ERROR);
		status = false;
	}

	if (params.endTime <= params.startTime || (params.startTime + params.timeStep) > params.endTime)
	{
		LogMan::Log("ERROR! End time must be greater than Start time with at least Time Step", LOG_ERROR);
		status = false;
	}

	//Check warn
	//TODO add warning here about low iteration time and such.


	return status;
}

bool LoadInputRasters(ModelParameters const & params)
{
	//Note Variable Precipitation rasters are not loaded initially, but loaded at runtime. A pass should, however, be carried initially
	//to create the time column of the time series vs raster path (so whe can tell when to load which raster). Caveats of this approach is
	//that, due to inefficient GeoTIFF parser, this would signifcantly lower the performance.
	//An alternative would be to preload a timeseries for each node (loop over rasters for each pair of coords). Caveats is that this would
	//need memory equivalent to at least Matrix_f32[nodeCount][tsLength] plus a double[rasterCount], where tsLength = raster count if all raster
	//fall within range of simulation time.
	//A third alternative solution is to create these per-node timeseries, but cache them to disk instead. Performance would be slower than memory
	//load, but faster than on-demand raster load. Memory would be much less than memory load.

	//DEM always loaded
	bool status;
	
	status = LoadRaster(params.demPath, &demID);

	if (params.variableManningCoefficients)
		status = LoadRaster(params.manningCoefficientRasterPath, &manningRasterID);
	
	if (!status) //error already logged with the function calls above.
	{
		UnloadRaster(demID);
		if (params.variableManningCoefficients)
			UnloadRaster(manningRasterID);
	}

	return status;
}

bool Simulate(ModelParameters const & params)
{
	if (!CheckParameters(params))
		return false;
	
	if (!LoadInputRasters(params))
		return false;

	ConstructGlobalConductanceMatrices(params.timeStep);
	ConstructGlobalCapacitanceMatrix();
	ConstructGlobalBetaMatrix();


	return true;
}

void TestSimulate(std::string const & nodesPath)
{
	std::cout << "\n test simulation start\n";
	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	LoadCoordinatePairsCSV(nodesPath, nodes);
	
	TestBoundingBox();
	Triangulate(nodes, &triangles, &boundaryNodes);
}
