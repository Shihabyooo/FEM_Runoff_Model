#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//TODO add a cleanup method to clear the allocated memory (superTriangles, rasters) when program closes.

std::unordered_map<int, Triangle> triangles;
Triangle superTriangles[2];
std::vector<Vector2D> nodes;
std::vector<int> boundaryNodes;
Vector2D nodesSW, nodesNE;
//TODO improve this
size_t exitNode = 0; //Assume first node to be exit node 

//Matrices and Vectors
Matrix_f64 globalC, globalPsiX, globalPsiY;
Vector_f64 nodeSlopeX, nodeSlopeY, nodeManning; //TODO consider pre-computing sqrt(slope-x)/manning and sqrt(slope-y)/manning and caching them instead
Vector_f64 heads;

Matrix_f64 const * dem = NULL, * manningRaster = NULL, * slopes = NULL, * fdr = NULL;

int demID, manningRasterID, slopesID, fdrID;

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

bool GenerateMesh(std::string const & nodesPath, double superTrianglePadding)
{
	LogMan::Log("Attempting to generate FEM mesh");

	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	if (!FileIO::LoadCoordinatePairsCSV(nodesPath, nodes))
		return false;

	ComputeBoundingBox();
	Triangulate(nodes, superTrianglePadding, &triangles, &boundaryNodes, superTriangles);

	LogMan::Log("Succesfully generated mesh of " + std::to_string(triangles.size()) + " triangles!", LOG_SUCCESS);
	return true;
}

Triangle const * GetElementContainingPoint(Vector2D const & pos)
{
	if (triangles.size() < 1)
		return NULL;

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
		if (it->second.ContainsPoint(pos))
			return &it->second;

	return NULL;
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

	/*if (!FileIO::FileExists(params.demPath))
	{
		LogMan::Log("ERROR! Must supply a terrain DEM for the region.", LOG_ERROR);
		status = false;
	}*/

	if (!FileIO::FileExists(params.slopesPath))
	{
		LogMan::Log("ERROR! Must supply a terrain Slopes for the region.", LOG_ERROR);
		status = false;
	}

	if (!FileIO::FileExists(params.fdrPath))
	{
		LogMan::Log("ERROR! Must supply a flow direction map (AGNPS format) for the region.", LOG_ERROR);
		status = false;
	}

	if (!params.variablePrecipitation && !params.unitTimeSeries.IsValid())
	{
		LogMan::Log("ERROR! Must provide a valid Time-series", LOG_ERROR);
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
	
	bool status = true;

	//status = status && FileIO::LoadRaster(params.demPath, &demID, dem);
	status = status && FileIO::LoadRaster(params.slopesPath, &slopesID, (void const **) &slopes);
	status = status && FileIO::LoadRaster(params.fdrPath, &fdrID, (void const **) &fdr);

	if (params.variableManningCoefficients)
		status = status && FileIO::LoadRaster(params.manningCoefficientRasterPath, &manningRasterID, (void const **) &manningRaster);

	if (!status) //error already logged with the function calls above.
	{
		FileIO::UnloadRaster(demID);
		FileIO::UnloadRaster(slopesID);
		if (params.variableManningCoefficients)
			FileIO::UnloadRaster(manningRasterID);
	}

	return status;
}

void UnloadAllRasters()
{
	if (dem != NULL)
		FileIO::UnloadRaster(demID);
	dem = NULL;

	if (manningRaster != NULL)
		FileIO::UnloadRaster(manningRasterID);
	manningRaster = NULL;

	if (slopes != NULL)
		FileIO::UnloadRaster(slopesID);
	slopes = NULL;

	if (fdr != NULL)
		FileIO::UnloadRaster(fdrID);
	fdr = NULL;
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

//returns negative valued VectorInt if error
std::pair<Vector2Int, double> SampleNearestPixel(Vector2D position, Matrix_f64 const * raster, int rasterID)
{
	//TODO  cache rastermapping params
	double ** tiePoints = NULL;
	double * pixelScale = NULL;
	bool isUTM;
	Vector2Int dimensions;
	int samples;

	if (!FileIO::GetRasterMappingParameters(rasterID, dimensions, samples, isUTM, &tiePoints, &pixelScale))
	{
		//Error already logged in GetRasterMappingParameters();
		return std::pair<Vector2Int, double>(Vector2Int(-1, -1), 0.0);
	}
	
	size_t _row, _column;
	double minDist = DBL_MAX;

	Vector2D ancrhoNE(tiePoints[1][0] + (pixelScale[0] / 2.0), tiePoints[1][1] - (pixelScale[1] / 2.0));
	for (size_t row = 0; row < raster->Rows(); row++)
		for (size_t column = 0; column < raster->Columns(); column++)
		{
			Vector2D pixelPos(	ancrhoNE.x + column * pixelScale[0],
								ancrhoNE.y - row * pixelScale[1]);

			double distanceToCentroid = pixelPos.DistanceTo(position);

			if (distanceToCentroid < minDist)
			{
				_row = row;
				_column = column;
				minDist = distanceToCentroid;
			}
		}

	if (isnan(raster->GetValue(_row, _column)))
	{
		LogMan::Log("Warning! Caught noValue pixel sampling for node ", LOG_WARN);
		return std::pair<Vector2Int, double>(Vector2Int(-1, -1), 0.0);
	}

	std::pair< Vector2Int, double > result(Vector2Int(_column, _row), raster->GetValue(_row, _column));

	//cleanup memory
	delete[] pixelScale;
	delete[] tiePoints[0];
	delete[] tiePoints[1];
	delete[] tiePoints;

	return result;
}

bool CacheManningCoefficients(ModelParameters const & params)
{
	nodeManning = Vector_f64(nodes.size());

	size_t counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		if (!params.variableManningCoefficients)
		{
			nodeManning[counter] = params.fixedManningCoeffient;
		}
		else
		{
			LogMan::Log("ERROR! Variable manning input not yet implemented.", LOG_ERROR);
			return false;
		}
		
		counter++;
	}

	return true;
}

bool CacheSlopes(ModelParameters const & params)
{
	//agnps fdr format: 1 = north, 2 = NE, 3 = East, 4 = SE, 5 = South, 6 = SW, 7 = West, 8 = NW
	
	//TODO revise this method so it creates a polygon connecting the centroids for each triangle that contains each tested vertex,\
	then sampling for pixels covered by this polygon and averaging the results (after factoring). Same for manning sampling.

	nodeSlopeX = Vector_f64(nodes.size());
	nodeSlopeY = Vector_f64(nodes.size());

	size_t counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		std::pair<Vector2Int, double> nodeSlopePixel = SampleNearestPixel(*it, slopes, slopesID);
		if (nodeSlopePixel.first.x < 0)
			LogMan::Log("Warning! No slopes sampled for node: " + std::to_string(counter), LOG_WARN);

		std::pair<Vector2Int, double> nodeFDRPixel = SampleNearestPixel(*it, fdr, fdrID);
		if (nodeFDRPixel.first.x < 0)
			LogMan::Log("Warning! No slopes sampled for node: " + std::to_string(counter), LOG_WARN);

		int dir = lround(nodeFDRPixel.second);

		if (dir == 1 || dir == 5)
			nodeSlopeY[counter] = nodeSlopePixel.second / 100.0;
		else if (dir == 3 || dir == 7)
			nodeSlopeX[counter] = nodeSlopePixel.second / 100.0;
		else
		{
			double component = nodeSlopePixel.second * 0.7071067811865475244 / 100.0; //adjacent = hypotenuse * cos(45)
			nodeSlopeX[counter] = component;
			nodeSlopeY[counter] = component;
		}

		counter++;
	}

	return true;
}

//Precipitation returned as meters per hour.
//current impl doesn't need triangle, but later it would.
double GetCurrentPrecipitation(double time, ModelParameters const & params, Triangle const & triangle)
{
	if (params.variablePrecipitation)
	{
		LogMan::Log("Warning! Variable precipitation is not yet implemented!", LOG_WARN);
		return 0.0;
	}
	else
	{
		return params.unitTimeSeries.SampleRate(time) / 1000.0;// , params.timeStep, params.precipitationTemporalInterpolationType);
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
		double newPrecipitation = GetCurrentPrecipitation(time + params.timeStep, params, it->second);
		double elementContrib = (it->second.area / 3.0) * (( 1.0 - params.femOmega) * it->second.elementPrecipitation + params.femOmega * newPrecipitation);
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
		
		for (int i = 0; i < 3; i++)
		{
			outVectorX[vert[i]] += sqrt(nodeSlopeX[vert[i]]) * pow(heads[vert[i]], 5.0 / 3.0) / nodeManning[vert[i]];
			outVectorY[vert[i]] += sqrt(nodeSlopeY[vert[i]]) * pow(heads[vert[i]], 5.0 / 3.0) / nodeManning[vert[i]];
		}
	}

	//test force q at bounderies to zero
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	{
		outVectorX[*it] = 0.0;
		outVectorY[*it] = 0.0;
	}
}

//test
Vector_f64 _oldHeads, _newHeads, _q_x_old, _q_x_new, _q_y_old, _q_y_new, _precipContrib, _RHS;
void TestShowInternalRHSVectors()
{
	Vector_f64 chold = globalC * _oldHeads;
	std::cout << "\n id | oldH  |  newH  ||| [C]{h0}|  qx_0  |  qx_1 |  qy_0  |  qy_1  | precip | RHS" << std::endl;
	for (size_t i = 0; i < _oldHeads.Rows(); i++)
		std::cout << std::fixed << std::setw(3) << i << " | " \
		<< std::setw(3) << std::setprecision(3) << _oldHeads[i] << " | " << std::setw(6) << _newHeads[i] << " ||| " \
		<< std::setw(6) << std::setprecision(1) << chold[i] << " | " \
		<< std::setw(6) << std::setprecision(3) << _q_x_old[i] << " | " << std::setw(3) << _q_x_new[i] << " | " \
		<< std::setw(6) << _q_y_old[i] << " | " << std::setw(6) << _q_y_new[i] << " | " \
		<< std::setw(6) << std::setprecision(1) << _precipContrib[i] << " | " << _RHS[i] << std::endl;
}
//endtest

void ComputeRHSVector(double time, ModelParameters const & params, Vector_f64 const & oldHeads, Vector_f64 const & newHeads, Vector_f64 & outRHS)
{
	//RHS is [C]{h_old} - dT [Psi-X] ((1-omega) * q_x_old + omega * q_x_new)
	//		- dT [Psi-Y] ((1-omega) * q_y_old + omega * q_y_new)
	//		+ PrecipitationVector

	Vector_f64 q_x_old, q_x_new, q_y_old, q_y_new;
	ComputeDischargeVectors(params, oldHeads, q_x_old, q_y_old);
	ComputeDischargeVectors(params, newHeads, q_x_new, q_y_new);

	//PsiX and PsiY already have dT multiplied with them.
	outRHS = globalC * oldHeads
			- ((globalPsiX) * ((q_x_old * (1 - params.femOmega)) + (q_x_new * params.femOmega)))
			- ((globalPsiY) * ((q_y_old * (1 - params.femOmega)) + (q_y_new * params.femOmega)))
			+ ComputePreciptationVector(time, params);
			
	//test
	//cache so we can display once after end of internal loop.
	_oldHeads = oldHeads;
	_newHeads = newHeads;
	_q_x_old = q_x_old;
	_q_x_new = q_x_new;
	_q_y_old = q_y_old;
	_q_y_new = q_y_new;
	_precipContrib = ComputePreciptationVector(time, params);
	_RHS = outRHS;
}

bool Simulate(ModelParameters const & params)
{
	LogMan::Log("Starting a simulation run");

	if (!CheckParameters(params))
		return false;
	
	UnloadAllRasters();
	if (!LoadInputRasters(params))
		return false;

	//Special consideration. Since the boundary node listing includes our exit node, we have to manually remove it.
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == exitNode)
		{
			boundaryNodes.erase(it);
			break;
		}

	LogMan::Log("Constructing global matrices and vectors");
	ConstructGlobalCapacitanceMatrix(params.useLumpedForm);
	ConstructGlobalPsiMatrices(params.timeStep);

	if (!CacheManningCoefficients(params))
		return false;
	
	if (!CacheSlopes(params))
		return false;

	//test
	std::cout << "\n===================================================\n";
	std::cout << "Global Capacitance";
	std::cout << "\n===================================================\n";
	globalC.DisplayOnCLI(0);

	std::cout << "\n===================================================\n";
	std::cout << "Global Psi-x";
	std::cout << "\n===================================================\n";
	globalPsiX.DisplayOnCLI(0);

	std::cout << "\n===================================================\n";
	std::cout << "Global Psi-y";
	std::cout << "\n===================================================\n";
	globalPsiY.DisplayOnCLI(0);
	
	std::cout << "\n===================================================\n";
	std::cout << "Boundery Nodes";
	std::cout << "\n===================================================\n";
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		std::cout << *it << std::endl;
	
	std::cout << "\n===================================================\n";
	std::cout << "Slopes and Manning roughness coef";
	std::cout << "\n===================================================\n";

	std::cout << "node |  n  |  Sx  |  Sy\n";
	for (size_t i = 0; i < nodeSlopeX.Rows(); i++)
		std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " << nodeManning[i] << " | " << nodeSlopeX[i] << " | " << nodeSlopeY[i] << std::endl;
		
	/*slopes->DisplayOnCLI();
	fdr->DisplayOnCLI();*/

	//return false;
	//end test

	//Set initial heads to zero (Dry conditions).
	heads = Vector_f64(nodes.size());

	//Capacitance matrix adjusted for boundary conditions
		//TODO do this when constructing this matrix.
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	{
		//https://finite-element.github.io/7_boundary_conditions.html
		for (size_t i = 0; i < globalC.Columns(); i++)
		{
			globalC[*it][i] = 0.0;
			/*globalPsiX[*it][i] = 0.0;
			globalPsiY[*it][i] = 0.0;*/
		}

		globalC[*it][*it] = 1.0;
		/*globalPsiX[*it][*it] = 1.0;
		globalPsiY[*it][*it] = 1.0;*/
	}

	//test
	/*for (size_t i = 0; i < globalC.Columns(); i++)
		globalC[exitNode][i] = 0.0;
	globalC[exitNode][exitNode] = 1.0;*/


	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop");
	
	while (time <= params.endTime)
	{
		std::cout << "\n\n===================================================\n";
		LogMan::Log("At T= " + std::to_string(time));
		std::cout << "===================================================\n";
		
		Vector_f64 newHeads = heads * 1.1;

		//internal loop
		for (size_t i = 0; i <= params.maxInternalIterations; i++)
		{
			Vector_f64 RHS;
			ComputeRHSVector(time, params, heads, newHeads, RHS);

			Vector_f64 residuals; //needed by solver
			Vector_f64 fixedNewH;

			//Adjust heads for boundary cond
			//https://finite-element.github.io/7_boundary_conditions.html
			for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
				RHS[*it] = 0.0; //this, plus the adjustment to globalC above, ensures resulting h for this node always = 0
			//RHS[exitNode] = 0.0;
			
			if (!Solve(globalC, RHS, fixedNewH, residuals, params))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}

			//force computed newHeads to be positive (not sure about this)
			for (size_t i = 0; i < fixedNewH.Rows(); i++)
				fixedNewH[i] = Max(fixedNewH[i], 0.0);

			if ((newHeads - fixedNewH).Magnitude() <= params.internalResidualTreshold)
			{
				newHeads = fixedNewH;
				break;
			}
			else if (i >= params.maxInternalIterations) //test
				std::cout << "reached maxInternalIterations without reaching appropriate h\n";
			
			newHeads = fixedNewH;
		}

		//test
		TestShowInternalRHSVectors();

		Vector_f64 qx, qy;
		ComputeDischargeVectors(params, heads, qx, qy);
		/*for (size_t i = 0; i < heads.Rows(); i++)
			std::cout << qx[i] << "\t" << qy[i] << std::endl;*/
		double area = 0.0;
		for (auto it = triangles.begin(); it != triangles.end(); ++it)
			if (it->second.ContainsVertex(exitNode))
			{
				area = it->second.area;
				break;
			}
		std::cout << "At time: " << time << "out discharge (x, y) : " << std::setprecision(4) << qx[exitNode] * sqrt(area) << ", " << qy[exitNode] * sqrt(area) << std::endl;
		//end test

		heads = newHeads;
		time += params.timeStep;
	}
	
	return true;
}