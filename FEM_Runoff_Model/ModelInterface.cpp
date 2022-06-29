#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//std::unordered_map<int, Triangle> triangles;
//std::vector<Vector2> nodes;
std::unordered_map<int, Triangle> triangles;
std::vector<Vector2D> nodes;
std::vector<int> boundaryNodes;
Vector2D nodesSW, nodesNE;

//Matrices and Vectors
Matrix_f32 globalC, globalPsiX, globalPsiY;
//Vector_f32 globalBeta;

Matrix_f32 const * dem = NULL, * manningRaster = NULL, * fdr;

int demID, manningRasterID, fdrID;

void ComputeBoundingBox()
{
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
	//[Psi-X_e] = 1/6 * delta T *	|	yj-yk	yk-yi	yi-yj	|
	//								|	yj-yk	yk-yi	yi-yj	|
	//								|	yj-yk	yk-yi	yi-yj	|

	//[Psi-Y_e] = 1/6 * delta T *	|	xk-xj	xi-xk	xj-xi	|
	//								|	xk-xj	xi-xk	xj-xi	|
	//								|	xk-xj	xi-xk	xj-xi	|
	//Where x(i/j/k) and y(i/j/k) are the x, y coord of the i/j/kth node.

	//Create global  Psi-X and Psi-y of size (node x nodes), zero initial value.
	//loop over each triangle in mesh
		//loop over permutations of vertIDs (3 verts per triangle = 9 permutations)
			//Psi-X[permutation] += Psi-X_e[localized permutation], e.g. if triangle is verts 3, 5, 8  (i, j, k)\
			Psi-x[3][3] += Psi-x_e[0][0] = 1/6 * dT * (yj-yk), and\
			Psi-x[5][8] += Psi-x_e[1][2] = 1/6 * dT * (yK-yI), etc
			//Ditto for Psi-Y
	
	globalPsiX = Matrix_f32(nodes.size(), nodes.size());
	globalPsiY = Matrix_f32(nodes.size(), nodes.size());

	double multiplier = deltaT / 6.0;

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		int const * vert = it->second.vertIDs; //to simplify lines bellow.

		//Since all values in a column (for each element matrix) have the same value, we compute them first hand then increment global matrix.
		double valX1 = multiplier * (it->second.nodes[1].y - it->second.nodes[2].y); //1/6 * deltaT * (yj - yk)
		double valX2 = multiplier * (it->second.nodes[2].y - it->second.nodes[0].y); //1/6 * deltaT * (yk - yi)
		double valX3 = multiplier * (it->second.nodes[0].y - it->second.nodes[1].y); //1/6 * deltaT * (yi - yj)

		double valY1 = multiplier * (it->second.nodes[2].x - it->second.nodes[1].x); //1/6 * deltaT * (xk - xj)
		double valY2 = multiplier * (it->second.nodes[0].x - it->second.nodes[2].x); //1/6 * deltaT * (xi - xk)
		double valY3 = multiplier * (it->second.nodes[2].x - it->second.nodes[0].x); //1/6 * deltaT * (xk - xi)

		//Increment global matrix at position define by node IDs (3^2 = 9 positions)
		globalPsiX[vert[0]][vert[0]] += valX1;
		globalPsiX[vert[0]][vert[1]] += valX2;
		globalPsiX[vert[0]][vert[2]] += valX3;

		globalPsiX[vert[1]][vert[0]] += valX1;
		globalPsiX[vert[1]][vert[1]] += valX2;
		globalPsiX[vert[1]][vert[2]] += valX3;

		globalPsiX[vert[2]][vert[0]] += valX1;
		globalPsiX[vert[2]][vert[1]] += valX2;
		globalPsiX[vert[2]][vert[2]] += valX3;
		
		globalPsiY[vert[0]][vert[0]] += valY1;
		globalPsiY[vert[0]][vert[1]] += valY2;
		globalPsiY[vert[0]][vert[2]] += valY3;

		globalPsiY[vert[1]][vert[0]] += valY1;
		globalPsiY[vert[1]][vert[1]] += valY2;
		globalPsiY[vert[1]][vert[2]] += valY3;

		globalPsiY[vert[2]][vert[0]] += valY1;
		globalPsiY[vert[2]][vert[1]] += valY2;
		globalPsiY[vert[2]][vert[2]] += valY3;
	}
}

void ConstructGlobalCapacitanceMatrix(bool isLumped)
{
	//Lumped Capacitance matrix for each element is 3x3 matrix
	//[C_e] = A/3 * |	1	0	0	|
	//				|	0	1	0	|
	//				|	0	0	1	|
	//Where A is the area of element.


	//Consistent Capacitance matrix for each element is 3x3 matrix
	//[C_e] = A/12 *	|	2	1	1	|
	//					|	1	2	1	|
	//					|	1	1	2	|
	//Where A is the area of element.

	globalC = Matrix_f32(nodes.size(), nodes.size());

	//since element matrix comprises of two distinct val, we compute them per element depending on whether lumped formulation or consistent.
	double diagVal, otherVal;

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		//compute diagVal and otherVal
		if (isLumped)
		{
			diagVal = it->second.area / 3.0; //(A/3.0) * 1.0
			otherVal = 0.0;
		}
		else
		{
			diagVal = 2.0 * it->second.area / 12.0; 
			otherVal = it->second.area / 12.0; //(A/12.0) * 1.0
		}
		
		int const * vert = it->second.vertIDs;

		globalC[vert[0]][vert[0]] += diagVal;
		globalC[vert[1]][vert[1]] += diagVal;
		globalC[vert[2]][vert[2]] += diagVal;

		globalC[vert[0]][vert[1]] += otherVal;
		globalC[vert[0]][vert[2]] += otherVal;

		globalC[vert[1]][vert[0]] += otherVal;
		globalC[vert[1]][vert[2]] += otherVal;

		globalC[vert[2]][vert[0]] += otherVal;
		globalC[vert[2]][vert[1]] += otherVal;
	}
}

//void ConstructGlobalBetaMatrix()
//{
//	//Beta matrix for each element is a 3x1 vector
//	//{Beta_e} = A/3 *	|	1	|
//	//					|	1	|
//	//					|	1	|
//	//Where A is the area of element.
//
//	////Construct global matrix similar to proceudre in ConstructGlobalConductanceMatrices(), difference being is that we're only working\
//	on a row level.
//
//	globalBeta = Vector_f32(nodes.size());
//
//	for (auto it = triangles.begin(); it != triangles.end(); ++it)
//	{
//		double areaThird = it->second.area / 3.0;
//		int const * vert = it->second.vertIDs; //to simplify lines bellow.
//
//		globalBeta[vert[0]] += areaThird;
//		globalBeta[vert[1]] += areaThird;
//		globalBeta[vert[2]] += areaThird;
//	}
//}

double SampleRegionAveragedValue(Triangle const & element, Matrix_f32 const * raster, int rasterID) //get the average value of pixels in raster covered by element
{
	//This is a rasterization/coverage test problem. But for now, we'll do it a very crude way, i.e. test if centre of pixel lies within triangle.
	//TODO improve this

	double ** tiePoints = NULL;
	double * pixelScale = NULL;
	bool isUTM;

	if (!FileIO::GetRasterMappingParameters(rasterID, isUTM, tiePoints, pixelScale))
	{
		//Error already logged in GetRasterMappingParameters();
		//LogMan::Log("ERROR! At loading raster mapping parameters", LOG_ERROR);

		return 0.0;
	}

	std::vector<double> pixelValues;

	for (size_t i = 0; i < raster->Rows(); i++)
		for (size_t j = 0; j < raster->Columns(); j++)
		{
			Vector2D pixelPos(	tiePoints[1][0] + j * pixelScale[0],
								tiePoints[1][1] - i * pixelScale[1]);

			if (element.ContainsPoint(pixelPos))
				pixelValues.push_back(raster->GetValue(i, j));
		}
	
	//There could be a case where element is small that it fits inside a pixel but not cover it's centre.
	if (pixelValues.size() < 1)
	{
		//In this case, simply compute the centroid (average of ndoes), and find pixel that contains centroid, return its value
		Vector2D centroid = element.Centroid();
		size_t row, column;
		double minDist = DBL_MAX;
		
		for (size_t i = 0; i < raster->Rows(); i++)
			for (size_t j = 0; j < raster->Columns(); j++)
			{
				Vector2D pixelPos(tiePoints[1][0] + j * pixelScale[0],
					tiePoints[1][1] - i * pixelScale[1]);

				double distanceToCentroid = pixelPos.DistanceTo(centroid);
				if (distanceToCentroid < minDist)
				{
					row = i;
					column = j;
				}
			}

		return raster->GetValue(row, column);
	}

	double sum = 0.0;
	for (auto it = pixelValues.begin(); it != pixelValues.end(); ++it)
		sum += *it;
	
	return sum / static_cast<double>(pixelValues.size());
}

void CacheManningCoefficients(ModelParameters const & params)
{
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		if (!params.variableManningCoefficients)
		{
			it->second.manningCoef = params.fixedManningCoeffient;
		}
		else
		{
			it->second.manningCoef = SampleRegionAveragedValue(it->second, manningRaster, manningRasterID);
		}
	}
}

//Matrix_f32 D8FDR() //computes Flow Direction Map in 8 direction. Encoded as follows:
//{
//
//}

void CacheSlopes(ModelParameters const & params)
{
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{

	}
}

bool CheckParameters(ModelParameters const & params)
{
	LogMan::Log("Checking parameters");
	bool status = true;

	//Check fails

	//TODO for raster file checks, check also that they are georeffed rasters and inclusive of the mesh boundary.
	if (triangles.size() < 1) //existence of tris implies existence of nodes (though I probably need to protect both from change outside this file)
	{
		LogMan::Log("ERROR! No loaded mesh.", LOG_ERROR);
		status = false;
	}

	if (!FileIO::FileExists(params.demPath))
	{
		LogMan::Log("ERROR! Must supply a terrain DEM for the region.", LOG_ERROR);
		status = false;
	}

	//if (params.variablePrecipitation && params.unitTimeSeries == NULL) //TODO check for time series rasters also goes here
	if (params.variablePrecipitation && !params.unitTimeSeries.IsValid()) //TODO check for time series rasters also goes here
	{
		LogMan::Log("ERROR! Invalid Time-series.", LOG_ERROR);
		status = false;
	}
	else if (!params.variablePrecipitation && params.fixedPrecipitationValue <= 0.0)
	{
		LogMan::Log("ERROR! Must set a positive precipitation value when using fixed precipitation", LOG_ERROR);
		status = false;
	}

	if (params.variableManningCoefficients && !FileIO::FileExists(params.manningCoefficientRasterPath))
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

	if (params.femOmega < 0.0 || params.femOmega > 1.0)
	{
		LogMan::Log("ERROR! Omega must be from 0.0 to 1.0", LOG_ERROR);
		status = false;
	}

	//Check warn
	//TODO add warning here about low iteration time and such.

	return status;
}

bool LoadInputRasters(ModelParameters const & params)
{
	LogMan::Log("Loading rasters");
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
	
	status = FileIO::LoadRaster(params.demPath, &demID, dem);

	if (params.variableManningCoefficients)
		status = FileIO::LoadRaster(params.manningCoefficientRasterPath, &manningRasterID, manningRaster);
	
	if (!status) //error already logged with the function calls above.
	{
		FileIO::UnloadRaster(demID);
		if (params.variableManningCoefficients)
			FileIO::UnloadRaster(manningRasterID);
	}

	return status;
}

bool GenerateMesh(std::string const & nodesPath)
{
	LogMan::Log("Attempting to generate FEM mesh");
	
	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	if (!FileIO::LoadCoordinatePairsCSV(nodesPath, nodes))
		return false;

	ComputeBoundingBox();
	Triangulate(nodes, &triangles, &boundaryNodes);

	//Validate Triangles
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
		if (!it->second.Validate())
			LogMan::Log("Warning! Invalid Triangle " + std::to_string(it->second.id) + ". Area: " + std::to_string(it->second.area), LOG_WARN);
		
	LogMan::Log("Succesfully generated mesh!", LOG_SUCCESS);
	return true;
}

bool LoadTimeSeries(std::string const & path, TimeSeries & ts)
{
	return FileIO::LoadTimeSeries(path, ts);
}

double GetCurrentPrecipitation(double time, ModelParameters const & params, Triangle const & triangle) //current impl doesn't need triangle, but later it would.
{
	if (params.variablePrecipitation)
	{
		return params.unitTimeSeries.Sample(time, params.precipitationTemporalInterpolationType);
	}
	else
	{
		return params.fixedPrecipitationValue;
	}
}

void ComputePreciptationVector(double time, ModelParameters const & params, Vector_f32 & outVector)
{
	//Preciptation Vector is the last term of the RHS of the formulation. i.e. dT * {Beta} * ((1 - omega) * Pe_t + omega * Pe_t+dt)
	
	//Beta matrix for each element is a 3x1 vector
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.


	//zero out outVector before doing anything.
	outVector = Vector_f32(nodes.size());

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		int const * vert = it->second.vertIDs; //to simplify lines bellow.
		double newPrecipitation = GetCurrentPrecipitation(time, params, it->second);
		double elementContrib = it->second.area / 3.0 * (( 1 - params.femOmega) * it->second.elementPrecipitation + params.femOmega * newPrecipitation);
		it->second.elementPrecipitation = newPrecipitation;

		outVector[vert[0]] += elementContrib;
		outVector[vert[1]] += elementContrib;
		outVector[vert[2]] += elementContrib;
	}
}

void ComputeDischargeVector(ModelParameters const & params, Vector_f32 const & heads, Vector_f32 & outVector)
{
	//Psi-X and Psi-Y (for each element) are 3x3 matrices.
	//[Psi-X_e] = 1/6 * delta T *	|	yj-yk	yk-yi	yi-yj	|
	//								|	yj-yk	yk-yi	yi-yj	|
	//								|	yj-yk	yk-yi	yi-yj	|

	//[Psi-Y_e] = 1/6 * delta T *	|	xk-xj	xi-xk	xj-xi	|
	//								|	xk-xj	xi-xk	xj-xi	|
	//								|	xk-xj	xi-xk	xj-xi	|
	//Where x(i/j/k) and y(i/j/k) are the x, y coord of the i/j/kth node.

	//For Qs, using Manning equation. Q = (A/n) * R^(2/3) * sqrt(S)
	//Wide channel assumption, R ~= h.

	//The per-Element qs vectors are given as:
	// dT * Psi-X * ((1-omega) q-X_t + omega * q-X_t+dT) + dT * Psi-Y * ((1-omega) q-Y_t + omega * q-Y_t+dT)
	
	
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		int const * vert = it->second.vertIDs; //to simplify lines bellow.
		double mult = (it->second.area / it->second.manningCoef);
		double multX = mult * sqrt(it->second.slopeX);
		double multY = mult * sqrt(it->second.slopeY);
		double newDischargeX[3] = { multX * pow(heads.GetValue(vert[0]), 2.0 / 3.0), multX * pow(heads.GetValue(vert[1]), 2.0 / 3.0), multX * pow(heads.GetValue(vert[2]), 2.0 / 3.0) };
		double newDischargeY[3] = { multY * pow(heads.GetValue(vert[0]), 2.0 / 3.0), multY * pow(heads.GetValue(vert[1]), 2.0 / 3.0), multY * pow(heads.GetValue(vert[2]), 2.0 / 3.0) };

	}
}

bool Simulate(ModelParameters const & params)
{
	LogMan::Log("Starting a simulation run");
	//TODO before simulating, unload all loaded rasters
	if (!CheckParameters(params))
		return false;
	
	if (!LoadInputRasters(params))
		return false;

	LogMan::Log("Constructing global matrices and vectors");
	ConstructGlobalConductanceMatrices(params.timeStep);
	ConstructGlobalCapacitanceMatrix(params.useLumpedForm);
	//ConstructGlobalBetaMatrix();
	CacheManningCoefficients(params);
	CacheSlopes(params);

	//Set initial heads to zero (Dry conditions).
	Vector_f32 head(nodes.size());

	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop at time: " + std::to_string(time));
	
	Vector_f32 prectipitationComponent(nodes.size());
	Vector_f32 dischargeComponent(nodes.size());

	while (time <= params.endTime)
	{
		LogMan::Log("At T= " + std::to_string(time));
		//Compute effective rainfall for each element
		ComputePreciptationVector(time, params, prectipitationComponent);

		//Compute qx and qy vectors using previous pass heads (for first loop, it's the initial head)
		ComputeDischargeVector(params, head, dischargeComponent);
		//Internal loop
			//Solve the system of equations.
			//Check the resulting heads with the ones assumed for qx and qy, if the difference is large, recompute qx and qy, loop again.
			//if difference is accepable, break internal loop.
		//handle data storage for results, residuals and any relative statistics, prepare for next loop.
		//Display results and return control to user.

		time += params.timeStep;
	}
	

	return true;
}