#pragma once
#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//TODO add a cleanup method to clear the allocated memory (superTriangles, rasters) when program closes.
std::unordered_map<int, Rectangle> rectangles;
std::unordered_map<int, Triangle> triangles;
Triangle superTriangles[2];
std::vector<Vector2D> nodes;
std::vector<int> boundaryNodes;
Vector2D nodesSW, nodesNE;
Vector2D shedSW, shedNE;
//TODO improve this
size_t exitNode = 0; //Assume first node to be exit node 
std::vector<Vector2D> shedBoundary;

//Matrices and Vectors
Matrix_f64 globalC, globalPsiX, globalPsiY;
Vector_f64 nodeSlopeX, nodeSlopeY, nodeManning; //TODO consider pre-computing sqrt(slope-x)/manning and sqrt(slope-y)/manning and caching them instead
Vector_f64 nodeFDR; //test
Vector_f64 nodeElevation;
Vector_f64 heads;

Matrix_f64 const * dem = NULL, * manningRaster = NULL, * slopes = NULL, * fdr = NULL;

int demID, manningRasterID, slopesID, fdrID;

void ComputeBoundingBox(std::vector<Vector2D> const & points, Vector2D & min, Vector2D & max)
{
	min = points[0];
	max = points[0];

	for (auto it = points.begin(); it < points.end(); it++)
	{
		Print(*it);
		min.x = Min(it->x, min.x);
		min.y = Min(it->y, min.y);
		max.x = Max(it->x, max.x);
		max.y = Max(it->y, max.y);
	}

	LogMan::Log("Loaded data with bounds: "
				+ std::to_string(min.x) + ", " + std::to_string(min.y) + " and "
				+ std::to_string(max.x) + ", " + std::to_string(max.y));
}

bool LoadWatershedBoundary(std::string const & boundaryPath)
{
	if (!FileIO::LoadVectorPath(boundaryPath, shedBoundary))
	{
		return false;
	}

	//test
	ComputeBoundingBox(shedBoundary, shedSW, shedNE);

	return true;
}

bool GenerateMesh(std::string const & nodesPath, double superTrianglePadding)
{
	LogMan::Log("Attempting to generate FEM mesh");

	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	if (!FileIO::LoadCoordinatePairsCSV(nodesPath, nodes))
		return false;

	ComputeBoundingBox(nodes, nodesSW, nodesNE);
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

bool UpdateNode(size_t id, Vector2D const & newPos)
{
	if (id >= nodes.size())
		return false;

	nodes[id] = newPos;

	//look for triangles using this node and update their areas, while logging warning if the change made it invalid
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		if (it->second.ContainsVertex(id))
		{
			for (int i = 0; i < 3; i++)
				it->second.nodes[i] = nodes[it->second.vertIDs[i]];

			if (!it->second.Validate())
				LogMan::Log("Warning! Updated triangle " + std::to_string(it->second.id) + " is invalid (Area: " + std::to_string(it->second.area) + ").", LOG_WARN);
		}
	}

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

	status = status && FileIO::LoadRaster(params.demPath, &demID, &dem);
	status = status && FileIO::LoadRaster(params.slopesPath, &slopesID, &slopes);
	status = status && FileIO::LoadRaster(params.fdrPath, &fdrID, &fdr);

	if (params.variableManningCoefficients)
		status = status && FileIO::LoadRaster(params.manningCoefficientRasterPath, &manningRasterID, &manningRaster);

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
	//test
	for (int i = 0; i < nodes.size(); i++)
	{
		std::string prefix = "internal";
		for (int j = 0; j < boundaryNodes.size(); j++)
			if (i == boundaryNodes[j])
				prefix = "external";
		
		std::cout << prefix.c_str() << std::fixed << " - " << i << " : " << std::setprecision(2) << globalPsiX[i][i] << " | " << std::setprecision(2) << globalPsiY[i][i] << std::endl;
	}
	
	//for (size_t i = 0; i < globalPsiX.Rows(); i++)
	//	for (size_t j = 0; j < globalPsiX.Columns(); j++)
	//	{
	//		globalPsiX[i][j] = abs(globalPsiX[i][j]);
	//		globalPsiY[i][j] = abs(globalPsiY[i][j]);
	//	}

	for (size_t i = 0; i < globalPsiX.Rows(); i++)
	{
		double rowSumX = 0.0, columnSumX = 0.0;
		double rowSumY = 0.0, columnSumY = 0.0;
		for (size_t j = 0; j < globalPsiX.Rows(); j++)
		{
			rowSumX += globalPsiX[i][j];
			columnSumX += globalPsiX[j][i];
			rowSumY += globalPsiY[i][j];
			columnSumY += globalPsiY[j][i];
		}
		std::cout << i << ":\t" << rowSumX << "\t" << columnSumX << "\t" << rowSumY << "\t" << columnSumY << "\n";
	}
	//endtest
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
	nodeFDR = Vector_f64(nodes.size());

	size_t counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		std::pair<Vector2Int, double> nodeSlopePixel = SampleNearestPixel(*it, slopes, slopesID);
		if (nodeSlopePixel.first.x < 0)
		{
			LogMan::Log("ERROR! No slopes sampled for node: " + std::to_string(counter), LOG_ERROR);
			return false;
		}

		std::pair<Vector2Int, double> nodeFDRPixel = SampleNearestPixel(*it, fdr, fdrID);
		if (nodeFDRPixel.first.x < 0)
		{
			LogMan::Log("ERROR! No slopes sampled for node: " + std::to_string(counter), LOG_ERROR);
			return false;
		}


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

		nodeFDR[counter] = dir;
		counter++;
	}

	return true;
}

bool CacheElevations()
{
	//TODO similar to CacheSlopes()

	nodeElevation = Vector_f64(nodes.size());
	size_t counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		std::pair<Vector2Int, double> nearestPixel = SampleNearestPixel(*it, dem, demID);
		if (nearestPixel.first.x < 0)
		{
			LogMan::Log("ERROR! No elevation sampled for node: " + std::to_string(counter), LOG_ERROR);
			return false;
		}
		nodeElevation[counter] = nearestPixel.second;
		counter++;
	}

	return true;
}

//void ZeroBoundaryNodes(Vector_f64 & targetVector)
//{
//	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
//		targetVector[*it] = 0.0;
//}

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

		//test
		if (!it->second.ContainsVertex(17) && !it->second.ContainsVertex(18) && !it->second.ContainsVertex(23) &&
			!it->second.ContainsVertex(24) && !it->second.ContainsVertex(28) && !it->second.ContainsVertex(29) &&
			!it->second.ContainsVertex(30) && !it->second.ContainsVertex(31) && !it->second.ContainsVertex(32) &&
			!it->second.ContainsVertex(33) )
		/*if (!it->second.ContainsVertex(47) && !it->second.ContainsVertex(48) && !it->second.ContainsVertex(49) &&
			!it->second.ContainsVertex(56) && !it->second.ContainsVertex(57) && !it->second.ContainsVertex(58) &&
			!it->second.ContainsVertex(59) && !it->second.ContainsVertex(60) && !it->second.ContainsVertex(68) &&
			!it->second.ContainsVertex(69) && !it->second.ContainsVertex(70) && !it->second.ContainsVertex(71) &&
			!it->second.ContainsVertex(77) && !it->second.ContainsVertex(78) && !it->second.ContainsVertex(79) && 
			!it->second.ContainsVertex(80) && !it->second.ContainsVertex(81) && !it->second.ContainsVertex(82) &&
			!it->second.ContainsVertex(83) && !it->second.ContainsVertex(84) && !it->second.ContainsVertex(85) &&
			!it->second.ContainsVertex(86))*/
			elementContrib = 0.0;
		//end test


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
			//double head = heads[vert[i]];
			double head = heads[vert[i]] - nodeElevation[vert[i]];

			double xContrib = sqrt(nodeSlopeX[vert[i]]) * pow(head, 5.0 / 3.0) / nodeManning[vert[i]];
			double yContrib = sqrt(nodeSlopeY[vert[i]]) * pow(head, 5.0 / 3.0) / nodeManning[vert[i]];
			
			/*double sign = head < 0 ? -1.0 : 1.0;
			double xContrib = sign * sqrt(nodeSlopeX[vert[i]]) * pow(abs(head), 5.0 / 3.0) / nodeManning[vert[i]];
			double yContrib = sign * sqrt(nodeSlopeY[vert[i]]) * pow(abs(head), 5.0 / 3.0) / nodeManning[vert[i]];*/

			outVectorX[vert[i]] += xContrib;
			outVectorY[vert[i]] += yContrib;

			//test
			/*int dir = lround(nodeFDR[vert[i]]);
			double signX = (dir == 2 || dir == 3 || dir == 4) ? 1.0 : -1.0;
			double signY = (dir == 1 || dir == 2 || dir == 8) ? 1.0 : -1.0;
			outVectorX[vert[i]] += signX * xContrib;
			outVectorY[vert[i]] += signY * yContrib;*/
			//end test

		}
	}
	/*for (auto i = 0; i < outVectorX.Rows(); i++)
		std::cout << outVectorX[i] << ", " << outVectorY[i] << std::endl;*/
	/*for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	{
		outVectorX[*it] = 0.0;
		outVectorY[*it] = 0.0;
	}*/
}

#pragma region Test
Vector_f64 _oldHeads, _newHeads, _q_x_old, _q_x_new, _q_y_old, _q_y_new, _precipContrib, _RHS;
std::vector<std::pair<double, double>> outletQs;
bool IsBoundaryNode(int id)
{
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == id)
			return true;	
	return false;
}

void TestShowInternalRHSVectors()
{
	Vector_f64 chold = globalC * _oldHeads;
	std::cout << "\n id   | oldH  |  newH  |||Diff||| [C]{h0}|  qx_0  |  qx_1 |  qy_0  |  qy_1  | precip | RHS" << std::endl;
	for (size_t i = 0; i < _oldHeads.Rows(); i++)
		std::cout << std::fixed << std::setw(3) << i \
		<< ( IsBoundaryNode(i) ? " B" : "  ") << " | "\
		<< std::setw(3) << std::setprecision(3) << _oldHeads[i] << " | " << std::setw(6) << _newHeads[i] << " ||| " \
		<< (_newHeads[i] > _oldHeads[i] ? "UP" : (_newHeads[i] == _oldHeads[i] ? "--" :  "DN")) << " ||| " \
		<< std::setw(6) << std::setprecision(1) << chold[i] << " | " \
		<< std::setw(6) << std::setprecision(3) << _q_x_old[i] << " | " << std::setw(3) << _q_x_new[i] << " | " \
		<< std::setw(6) << _q_y_old[i] << " | " << std::setw(6) << _q_y_new[i] << " | " \
		<< std::setw(6) << std::setprecision(1) << _precipContrib[i] << " | " << _RHS[i] << std::endl;
}
#pragma endregion

#pragma region Decomposed Approach
bool IsLowestVert(int localVertID, Triangle const & element)
{
	for (int i = 0; i < 3; i++)
		if (nodeElevation[element.vertIDs[localVertID]] > nodeElevation[element.vertIDs[i]])
			return false;
	return true;
}

void ComputeUVComponents(ModelParameters const & params, Vector_f64 const & oldHeads, Vector_f64 const & newHeads, Vector_f64 & outVectorX, Vector_f64 & outVectorY)
{
	//this approach breaks q_x (in dq_x/dx) to u*h, where u is velocity in x direction. The differential is broken to two term,  u*dh/dx + h*du/dx
	//the approx solution (for terms relating to q_x) becomes [u][Psi_x]{h} + [h][Psi_X]{u}
	//											|	ui	0	0	|
	//Where [u] is a diagonal matrix of form =	|	0	uj	0	|
	//											|	0	0	uk	|
	//Same for [h].
	//Psi matrix is the unchanged.
	// u itself is computed using manning, u = sqrt(Sx) * h^(2/3) / n

	//outVectorX = Vector_f64(nodes.size());
	//outVectorY = Vector_f64(nodes.size());
	
	Matrix_f64 uMat(nodes.size(), nodes.size());
	Matrix_f64 vMat(nodes.size(), nodes.size());
	Matrix_f64 hMat(nodes.size(), nodes.size());
	Vector_f64 uVec(nodes.size());
	Vector_f64 vVec(nodes.size());
	Vector_f64 hVec(nodes.size());

	Matrix_f64 absPsiX = globalPsiX;//test
	Matrix_f64 absPsiY = globalPsiY;//test
	/*for (int i = 0; i < nodes.size(); i++)
	{
		for (int j = 0; j < nodes.size(); j++)
		{
			absPsiX[i][j] = abs(absPsiX[i][j]);
			absPsiY[i][j] = abs(absPsiY[i][j]);
		}
	}*/

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		int const * verts = it->second.vertIDs;

		for (int i = 0; i < 3; i++)
		{
			int vert = verts[i];
			double alpha = sqrt(nodeSlopeX[verts[i]]) / nodeManning[verts[i]];
			double beta = sqrt(nodeSlopeY[verts[i]]) / nodeManning[verts[i]];
			/*bool isLowest = IsLowestVert(i, it->second);
			alpha = isLowest ? alpha : -1.0 * alpha;
			beta = isLowest ? beta : -1.0 * beta;*/

			double centralHeadDiff = (1.0 - params.femOmega) * (oldHeads[vert] - nodeElevation[vert]) + params.femOmega * (newHeads[vert] - nodeElevation[vert]);
			double centralHeadPoweredDiff = (1.0 - params.femOmega) * pow(oldHeads[vert] - nodeElevation[vert], 2.0 / 3.0) + params.femOmega * pow(newHeads[vert] - nodeElevation[vert], 2.0 / 3.0);
			uMat[vert][vert] += alpha * centralHeadPoweredDiff;
			vMat[vert][vert] += beta * centralHeadPoweredDiff;
			hMat[vert][vert] += centralHeadDiff;
			uVec[vert] += alpha * centralHeadPoweredDiff;
			vVec[vert] += beta * centralHeadPoweredDiff;
			hVec[vert] += centralHeadDiff;
		}
	}

	for (int i = 0; i < nodes.size(); i++)
		std::cout << i << " : " << hVec[i] << " | " << uVec[i] << " | " << hVec[i] << std::endl;

	/*outVectorX = (uMat * globalPsiX * hVec) + (hMat * globalPsiX * uVec);
	outVectorY = (vMat * globalPsiY * hVec) + (hMat * globalPsiY * vVec);*/

	outVectorX = (uMat * absPsiX * hVec) + (hMat * absPsiX * uVec);
	outVectorY = (vMat * absPsiY * hVec) + (hMat * absPsiY * vVec);
	
	std::cout << "---=-=-=-=-=\n";
	//outVectorX.DisplayOnCLI();
	//outVectorY.DisplayOnCLI();
}
#pragma endregion

Vector_f64 q_x_old, q_y_old;
Vector_f64 last_q_x_new, last_q_y_new;

void ComputeRHSVector(double time, ModelParameters const & params, Vector_f64 const & oldHeads, Vector_f64 const & newHeads, Vector_f64 & outRHS)
{
	//RHS is [C]{h_old} - dT [Psi-X] ((1-omega) * q_x_old + omega * q_x_new)
	//		- dT [Psi-Y] ((1-omega) * q_y_old + omega * q_y_new)
	//		+ PrecipitationVector

	Vector_f64 q_x_new, q_y_new;
	ComputeDischargeVectors(params, newHeads, q_x_new, q_y_new);

	//PsiX and PsiY already have dT multiplied with them.
	//TODO [C]{h0} and {P} don't change from (internal) iteration to the next. Should refactor this function to have one called\
	every time loop (external loop), and the other every internal loop, which add results of external loop to {q} vectors.
	/*outRHS = globalC * oldHeads
			- ((globalPsiX) * ((q_x_old * (1.0 - params.femOmega)) + (q_x_new * params.femOmega)))
			- ((globalPsiY) * ((q_y_old * (1.0 - params.femOmega)) + (q_y_new * params.femOmega)))
			+ ComputePreciptationVector(time, params);*/

	//Using UV decomp
	Vector_f64 uComp, vComp;
	ComputeUVComponents(params, oldHeads, newHeads, uComp, vComp);
	outRHS = (globalC * oldHeads)
		- uComp
		- vComp
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

	//Vector_f64 compX1 = ((globalPsiX) * ((q_x_old * (1.0 - params.femOmega)) + (q_x_new * params.femOmega)));
	//Vector_f64 compY1 = ((globalPsiY) * ((q_y_old * (1.0 - params.femOmega)) + (q_y_new * params.femOmega)));
	//for (int i = 0; i < nodes.size(); i++)
	//	std::cout << q_x_new[i] << " - " << q_y_new[i] << " | " << compX1[i] << " - " << compY1[i] << std::endl;
	//end test

	//TestShowInternalRHSVectors();
	//cache the newly computed qs for use in next step.
	last_q_x_new = std::move(q_x_new);
	last_q_y_new = std::move(q_y_new);
}

bool Simulate(ModelParameters const & params)
{
	LogMan::Log("Starting a simulation run");

#pragma region Init
	if (!CheckParameters(params))
		return false;
	
	UnloadAllRasters();

	if (!LoadInputRasters(params))
		return false;

	//boundaryNodes.push_back(39); //test

	//Special consideration. Since the boundary node listing includes our exit node, we have to manually remove it.
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == exitNode)
		{
			boundaryNodes.erase(it);
			break;
		}

	if (!CacheManningCoefficients(params))
		return false;

	if (!CacheSlopes(params))
		return false;

	if (!CacheElevations())
		return false;

	LogMan::Log("Constructing global matrices and vectors");
	ConstructGlobalCapacitanceMatrix(params.useLumpedForm);
	ConstructGlobalPsiMatrices(params.timeStep);

#pragma endregion

#pragma region Test
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
	std::cout << "Time series";
	std::cout << "\n===================================================\n";
	for (int i = 0; i < params.unitTimeSeries.size; i++)
		std::cout << params.unitTimeSeries.series[i].first << " - " << params.unitTimeSeries.series[i].second << std::endl;

	std::cout << "\n===================================================\n";
	std::cout << "Boundery Nodes";
	std::cout << "\n===================================================\n";
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		std::cout << *it << std::endl;
	
	std::cout << "\n===================================================\n";
	std::cout << "Elevations, Slopes, Manning roughness coef";
	std::cout << "\n===================================================\n";

	std::cout << "node | Elev  |   n  |  Sx  |  Sy\n";
	for (size_t i = 0; i < nodeSlopeX.Rows(); i++)
		std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " << nodeElevation[i] << " | "  << nodeManning[i] << " | " << nodeSlopeX[i] << " | " << nodeSlopeY[i] << std::endl;
	//return false;
#pragma endregion

	//Set initial heads to zero (Dry conditions).
	//heads = Vector_f64(nodes.size());
	
	//nodeElevation = Vector_f64(nodes.size()); //test .
	heads = nodeElevation;
	//init cached old qs (not directly. value of last_q_x_new will be moved to q_x_old at begining of every time step.
	last_q_x_new = Vector_f64(nodes.size());
	last_q_y_new = Vector_f64(nodes.size());

	//Capacitance matrix adjusted for boundary conditions
	Matrix_f64 adjustedC = globalC;
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	{
		//https://finite-element.github.io/7_boundary_conditions.html
		for (size_t i = 0; i < globalC.Columns(); i++)
			adjustedC[*it][i] = 0.0;

		adjustedC[*it][*it] = 1.0;
	}

	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop");
	
	while (time <= params.endTime)
	{
		std::cout << "\n\n===================================================\n";
		LogMan::Log("At T= " + std::to_string(time));
		std::cout << "===================================================\n";
	
		Vector_f64 newHeads = heads;// +Vector_f64(nodes.size(), 0.05);

		q_x_old = std::move(last_q_x_new);
		q_y_old = std::move(last_q_y_new);

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
				RHS[*it] = nodeElevation[*it]; //this, plus the adjustment to globalC above, ensures resulting h for this node always = 0

			//Test
			//{
			//	//double area = 0.0;
			//	double widthX = 0.0, widthY = 0.0;
			//	for (auto it = triangles.begin(); it != triangles.end(); ++it)
			//		if (it->second.ContainsVertex(exitNode))
			//		{
			//			//area = it->second.area;
			//			int const * verts = it->second.vertIDs;
			//			int otherVerts[2];
			//			otherVerts[0] = verts[0] == exitNode ? verts[1] : verts[0];
			//			otherVerts[1] = verts[1] == exitNode ? verts[2] : (verts[1] == otherVerts[0] ? verts[2] : verts[1]);
			//			
			//			widthX = abs(nodes[otherVerts[0]].y - nodes[otherVerts[1]].y);
			//			widthY = abs(nodes[otherVerts[0]].x - nodes[otherVerts[1]].x);
			//			break;
			//		}
			//	//RHS[exitNode] -= ((q_x_old[exitNode] + last_q_x_new[exitNode]) / 2.0 + (q_y_old[exitNode] + last_q_y_new[exitNode])/ 2.0) * sqrt(area) * params.timeStep;
			//	double Qx = q_x_old[exitNode] * widthX * params.timeStep / 6.0;
			//	double Qy = q_y_old[exitNode] * widthY * params.timeStep / 6.0;
			//	RHS[exitNode] -= (Qx + Qy);
			//}
			//end test

			/*if (!Solve(adjustedC, RHS, fixedNewH, residuals, params))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}*/
			//test
			fixedNewH = Vector_f64(nodes.size());
			for (size_t j = 0; j < RHS.Rows(); j++)
				fixedNewH[j] = RHS[j] / adjustedC[j][j];
			//end test
			
			//(fixedNewH - nodeElevation).DisplayOnCLI(10);
			//force computed newHeads to be positive (not sure about this)
			for (size_t j = 0; j < fixedNewH.Rows(); j++)
				fixedNewH[j] = Max(fixedNewH[j], nodeElevation[j]);
			

			if ((newHeads - fixedNewH).Magnitude() <= params.internalResidualTreshold)
			{
				newHeads = fixedNewH;
				break;
			}
			else if (i >= params.maxInternalIterations)
				LogMan::Log("Reached maxInternalIterations without reaching appropriate h", LOG_WARN);
			
			newHeads = fixedNewH;
		}

		//test
		TestShowInternalRHSVectors();
		double area = 0.0;
		double widthX = 0.0, widthY = 0.0;
		for (auto it = triangles.begin(); it != triangles.end(); ++it)
			if (it->second.ContainsVertex(exitNode))
			{
				//area = it->second.area;
				int const * verts = it->second.vertIDs;
				int otherVerts[2];
				otherVerts[0] = verts[0] == exitNode ? verts[1] : verts[0];
				otherVerts[1] = verts[1] == exitNode ? verts[2] : (verts[1] == otherVerts[0] ? verts[2] : verts[1]);

				widthX = abs(nodes[otherVerts[0]].y - nodes[otherVerts[1]].y);
				widthY = abs(nodes[otherVerts[0]].x - nodes[otherVerts[1]].x);
				break;
			}
		//RHS[exitNode] -= ((q_x_old[exitNode] + last_q_x_new[exitNode]) / 2.0 + (q_y_old[exitNode] + last_q_y_new[exitNode])/ 2.0) * sqrt(area) * params.timeStep;
		double Qx = q_x_old[exitNode] * widthX;
		double Qy = q_y_old[exitNode] * widthY;
		//std::cout << "At time: " << time << " out discharge (x, y) : " << std::setprecision(4) << q_x_old[exitNode] * sqrt(area) << ", " << q_y_old[exitNode] * sqrt(area) << std::endl;
		std::cout << "At time: " << time << " out discharge (x, y) : " << std::setprecision(4) << Qx << ", " << Qy << std::endl;
		outletQs.push_back(std::pair(Qx, Qy));
		//end test

		heads = newHeads;
		time += params.timeStep;

		//std::cout << "\nTestEnd";//test
		//return false;
	}

	std::cout << "Results:\ntime  |  Q_x   |   Q_y  |  Product\n";
	double t = 0.0;
	for (auto it = outletQs.begin(); it != outletQs.end(); ++it)
	{
		std::cout << std::fixed << std::setw(4) << t << " | "  <<  std::setprecision(2) << it->first << "  |  " << std::setprecision(2) << it->second 
			<< "  |  " << sqrt(it->first * it->first + it->second * it->second) << std::endl;
		t += params.timeStep;
	}

	return true;
}
