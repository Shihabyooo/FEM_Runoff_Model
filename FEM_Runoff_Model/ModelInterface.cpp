#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//std::unordered_map<int, Triangle> triangles;
//std::vector<Vector2> nodes;
std::unordered_map<int, Triangle> triangles;
std::vector<Vector2D> nodes;
std::vector<int> boundaryNodes;
Vector2D nodesSW, nodesNE;
size_t exitNode;

//Matrices and Vectors
Matrix_f64 globalC, globalPsiX, globalPsiY;
//Vector_f64 globalBeta;
Vector_f64 heads;

Matrix_f64 const * dem = NULL, * manningRaster = NULL, * fdr;

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

//Model-specific functions

//void ConstructGlobalConductanceMatrices(double deltaT)
//{
//	//Psi-X and Psi-Y (for each element) are 3x3 matrices.
//	//[Psi-X_e] = 1/6 * delta T *	|	yj-yk	yk-yi	yi-yj	|
//	//								|	yj-yk	yk-yi	yi-yj	|
//	//								|	yj-yk	yk-yi	yi-yj	|
//
//	//[Psi-Y_e] = 1/6 * delta T *	|	xk-xj	xi-xk	xj-xi	|
//	//								|	xk-xj	xi-xk	xj-xi	|
//	//								|	xk-xj	xi-xk	xj-xi	|
//	//Where x(i/j/k) and y(i/j/k) are the x, y coord of the i/j/kth node.
//
//	//Create global  Psi-X and Psi-y of size (node x nodes), zero initial value.
//	//loop over each triangle in mesh
//		//loop over permutations of vertIDs (3 verts per triangle = 9 permutations)
//			//Psi-X[permutation] += Psi-X_e[localized permutation], e.g. if triangle is verts 3, 5, 8  (i, j, k)\
//			Psi-x[3][3] += Psi-x_e[0][0] = 1/6 * dT * (yj-yk), and\
//			Psi-x[5][8] += Psi-x_e[1][2] = 1/6 * dT * (yK-yI), etc
//			//Ditto for Psi-Y
//	
//	globalPsiX = Matrix_f64(nodes.size(), nodes.size());
//	globalPsiY = Matrix_f64(nodes.size(), nodes.size());
//
//	double multiplier = deltaT / 6.0;
//
//	for (auto it = triangles.begin(); it != triangles.end(); ++it)
//	{
//		int const * vert = it->second.vertIDs; //to simplify lines bellow.
//
//		//Since all values in a column (for each element matrix) have the same value, we compute them first hand then increment global matrix.
//		double valX1 = multiplier * (it->second.nodes[1].y - it->second.nodes[2].y); //1/6 * deltaT * (yj - yk)
//		double valX2 = multiplier * (it->second.nodes[2].y - it->second.nodes[0].y); //1/6 * deltaT * (yk - yi)
//		double valX3 = multiplier * (it->second.nodes[0].y - it->second.nodes[1].y); //1/6 * deltaT * (yi - yj)
//
//		double valY1 = multiplier * (it->second.nodes[2].x - it->second.nodes[1].x); //1/6 * deltaT * (xk - xj)
//		double valY2 = multiplier * (it->second.nodes[0].x - it->second.nodes[2].x); //1/6 * deltaT * (xi - xk)
//		double valY3 = multiplier * (it->second.nodes[2].x - it->second.nodes[0].x); //1/6 * deltaT * (xk - xi)
//
//		//Increment global matrix at position define by node IDs (3^2 = 9 positions)
//		globalPsiX[vert[0]][vert[0]] += valX1;
//		globalPsiX[vert[0]][vert[1]] += valX2;
//		globalPsiX[vert[0]][vert[2]] += valX3;
//
//		globalPsiX[vert[1]][vert[0]] += valX1;
//		globalPsiX[vert[1]][vert[1]] += valX2;
//		globalPsiX[vert[1]][vert[2]] += valX3;
//
//		globalPsiX[vert[2]][vert[0]] += valX1;
//		globalPsiX[vert[2]][vert[1]] += valX2;
//		globalPsiX[vert[2]][vert[2]] += valX3;
//		
//		globalPsiY[vert[0]][vert[0]] += valY1;
//		globalPsiY[vert[0]][vert[1]] += valY2;
//		globalPsiY[vert[0]][vert[2]] += valY3;
//
//		globalPsiY[vert[1]][vert[0]] += valY1;
//		globalPsiY[vert[1]][vert[1]] += valY2;
//		globalPsiY[vert[1]][vert[2]] += valY3;
//
//		globalPsiY[vert[2]][vert[0]] += valY1;
//		globalPsiY[vert[2]][vert[1]] += valY2;
//		globalPsiY[vert[2]][vert[2]] += valY3;
//	}
//}

bool LoadInputRasters(ModelParameters const & params)
{
	LogMan::Log("Loading rasters");
	//Note Variable Precipitation rasters are not loaded initially, but loaded at runtime. A pass should, however, be carried initially
	//to create the time column of the time series vs raster path (so whe can tell when to load which raster). Caveats of this approach is
	//that, due to inefficient GeoTIFF parser, this would signifcantly lower the performance.
	//An alternative would be to preload a timeseries for each node (loop over rasters for each pair of coords). Caveats is that this would
	//need memory equivalent to at least Matrix_f64[nodeCount][tsLength] plus a double[rasterCount], where tsLength = raster count if all raster
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

	globalC = Matrix_f64(nodes.size(), nodes.size());

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

void ConstructGlobalPsiMatrices(double timeStep)
{
	//Psi-X and Psi-Y (for each element) are 3x3 matrices.
	//[Psi-X_e] = 1/6 * delta T *	|	yj-yk	yk-yi	yi-yj	|
	//								|	yj-yk	yk-yi	yi-yj	|
	//								|	yj-yk	yk-yi	yi-yj	|

	//[Psi-Y_e] = 1/6 * delta T *	|	xk-xj	xi-xk	xj-xi	|
	//								|	xk-xj	xi-xk	xj-xi	|
	//								|	xk-xj	xi-xk	xj-xi	|
	//Where x(i/j/k) and y(i/j/k) are the x, y coord of the i/j/kth node.


	//Construct Global Psi-x and Psi-y

	globalPsiX = Matrix_f64(nodes.size(), nodes.size());
	globalPsiY = Matrix_f64(nodes.size(), nodes.size());

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		const int * verts = it->second.vertIDs;
		const Vector2D & i = it->second.nodes[0];
		const Vector2D & j = it->second.nodes[1];
		const Vector2D & k = it->second.nodes[2];

		double contribX1 = (timeStep / 6.0) * (j.y - k.y);
		double contribX2 = (timeStep / 6.0) * (k.y - i.y);
		double contribX3 = (timeStep / 6.0) * (i.y - j.y);

		double contribY1 = (timeStep / 6.0) * (k.x - j.x);
		double contribY2 = (timeStep / 6.0) * (i.x - k.x);
		double contribY3 = (timeStep / 6.0) * (j.x - i.x);

		for (int i = 0; i < 3; i++)
		{
			globalPsiX[verts[i]][verts[0]] += contribX1;
			globalPsiX[verts[i]][verts[1]] += contribX2;
			globalPsiX[verts[i]][verts[2]] += contribX3;

			globalPsiY[verts[i]][verts[0]] += contribY1;
			globalPsiY[verts[i]][verts[1]] += contribY2;
			globalPsiY[verts[i]][verts[2]] += contribY3;
		}
	}
}

double SampleRegionAveragedValue(Triangle const & element, Matrix_f64 const * raster, int rasterID) //get the average value of pixels in raster covered by element
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

void CacheSlopes(ModelParameters const & params)
{
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{

	}
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

Vector_f64 ComputePreciptationVector(double time, ModelParameters const & params)
{
	//Preciptation Vector is the last term of the RHS of the formulation. i.e. dT * {Beta} * ((1 - omega) * Pe_t + omega * Pe_t+dt)
	
	//Beta matrix for each element is a 3x1 vector
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.

	//zero out outVector before doing anything.
	Vector_f64 result(nodes.size());

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		int const * vert = it->second.vertIDs; //to simplify lines bellow.
		double newPrecipitation = GetCurrentPrecipitation(time, params, it->second);
		double elementContrib = it->second.area / 3.0 * (( 1 - params.femOmega) * it->second.elementPrecipitation + params.femOmega * newPrecipitation);
		elementContrib *= params.timeStep;
		it->second.elementPrecipitation = newPrecipitation;

		result[vert[0]] += elementContrib;
		result[vert[1]] += elementContrib;
		result[vert[2]] += elementContrib;
	}
	return result;
}

void ComputeDischargeVectors(ModelParameters const & params, Vector_f64 const & heads, Vector_f64 & outVectorX, Vector_f64 & outVectorY)
{
	//For Qs, using Manning equation. Q = (A/n) * R^(2/3) * sqrt(S)
	//Wide channel assumption, R ~= h. And for q (Q per unit width) A = h
	//q = (sqrt(s) / n) * h ^ (5/3)
	
	outVectorX = Vector_f64(nodes.size());
	outVectorY = Vector_f64(nodes.size());
	
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		int const * vert = it->second.vertIDs; //to simplify lines bellow.
		double multX = sqrt(it->second.slopeX) / it->second.manningCoef;
		double multY = sqrt(it->second.slopeY) / it->second.manningCoef;

		for (int i = 0; i < 3; i++)
		{
			outVectorX[vert[i]] += multX * pow(heads[vert[i]], 5.0 / 3.0);
			outVectorY[vert[i]] += multY * pow(heads[vert[i]], 5.0 / 3.0);
		}
	}
}

void ComputeRHSVector(double time, ModelParameters const & params, Vector_f64 const & oldHeads, Vector_f64 const & newHeads, Vector_f64 & outRHS)
{
	//RHS is [C]{h_old} - dT [Psi-X] ((1-omega) * q_x_old + omega * q_x_new)
	//		- dT [Psi-Y] ((1-omega) * q_y_old + omega * q_y_new)
	//		+ PrecipitationVector

	Vector_f64 q_x_old, q_x_new, q_y_old, q_y_new;
	ComputeDischargeVectors(params, oldHeads, q_x_old, q_y_old);
	ComputeDischargeVectors(params, newHeads, q_x_new, q_y_new);

	outRHS = globalC * oldHeads
			- ((globalPsiX * params.timeStep) * ((q_x_old * (1 - params.femOmega)) + (q_x_new * params.femOmega)))
			- ((globalPsiY * params.timeStep) * ((q_y_old * (1 - params.femOmega)) + (q_y_new * params.femOmega)))
			+ ComputePreciptationVector(time, params);
			
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
	ConstructGlobalCapacitanceMatrix(params.useLumpedForm);
	ConstructGlobalPsiMatrices(params.timeStep);

	CacheManningCoefficients(params);
	CacheSlopes(params);

	//Set initial heads to zero (Dry conditions).
	heads = Vector_f64(nodes.size());

	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop at time: " + std::to_string(time));
	
	/*Vector_f64 prectipitationComponent(nodes.size());
	Vector_f64 dischargeComponent(nodes.size());*/

	while (time <= params.endTime)
	{
		LogMan::Log("At T= " + std::to_string(time));

		Vector_f64 newHeads = heads * 1.1;

		//Capacitance matrix adjusted for boundary conditions
		Matrix_f64 adjustedGlobalC = globalC;
		for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		{
			//https://finite-element.github.io/7_boundary_conditions.html
			for (size_t i = 0; i < adjustedGlobalC.Columns(); i++)
				adjustedGlobalC[*it][i] = 0.0;

			adjustedGlobalC[*it][*it] = 1.0;
		}

		//internal loop
		for (size_t i = 0; i <= params.maxInternalIterations; i++)
		{
			Vector_f64 RHS;
			ComputeRHSVector(time, params, heads, newHeads, RHS);

			//Adjust RHS for boundary cond
			for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
			{
				//https://finite-element.github.io/7_boundary_conditions.html
				RHS[*it] = 0.0;
			}

			Vector_f64 residuals; //needed by solver
			Vector_f64 fixedNewH;

			Solve(adjustedGlobalC, RHS, fixedNewH, residuals, params);

			if ((newHeads - fixedNewH).Magnitude() <= params.internalResidualTreshold)
			{
				std::cout << "reached appropriate h\n"; //test
				break;
			}
			else if (i >= params.maxInternalIterations) //test
				std::cout << "reached maxInternalIterations without reaching appropriate h\n";

			newHeads = fixedNewH;
		}

		//test
		Vector_f64 qx, qy;
		ComputeDischargeVectors(params, heads, qx, qy);
		std::cout << "At time: " << time << "out discharge (x, y) : " << qx[exitNode] << ", " << qy[exitNode] << std::endl;
		//end test

		heads = newHeads;
		time += params.timeStep;
	}
	
	return true;
}