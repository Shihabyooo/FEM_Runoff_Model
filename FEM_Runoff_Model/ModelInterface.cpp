#pragma once
#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//TODO add a cleanup method to clear the allocated memory (superTriangles, rasters) when program closes.
std::unordered_map<size_t, Rectangle> rectangles;
std::unordered_map<size_t, Triangle> triangles;
Vector2D superTriangles[6];
std::vector<Vector2D> nodes;
std::vector<size_t> boundaryNodes;
Vector2D nodesSW, nodesNE;
Vector2D shedSW, shedNE;
std::vector<Vector2D> shedBoundary;

//Matrices and Vectors
Matrix_f64 globalC;
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
		//Print(*it);
		min.x = Min(it->x, min.x);
		min.y = Min(it->y, min.y);
		max.x = Max(it->x, max.x);
		max.y = Max(it->y, max.y);
	}

	LogMan::Log("Loaded data with bounds: "
				+ std::to_string(min.x) + ", " + std::to_string(min.y) + " and "
				+ std::to_string(max.x) + ", " + std::to_string(max.y));
}

void ResetMeshs()
{
	nodes.clear();
	triangles.clear();
	rectangles.clear();
	boundaryNodes.clear();
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

	ResetMeshs();

	if (!FileIO::LoadCoordinatePairsCSV(nodesPath, nodes))
		return false;

	ComputeBoundingBox(nodes, nodesSW, nodesNE);
	if (!Triangulate(nodes, superTrianglePadding, &triangles, &boundaryNodes, superTriangles))
	{
		LogMan::Log("Failed to generate FEM mesh", LOG_ERROR);
		ResetMeshs();
		return false;
	}

	LogMan::Log("Succesfully generated mesh of " + std::to_string(triangles.size()) + " triangles!", LOG_SUCCESS);
	return true;
}

bool GenerateGridMesh(size_t resolution, double internalPadding, double raycastPadding)
{
	LogMan::Log("Attempting to generate FEM mesh");
	
	ResetMeshs();

	if (!GenerateGrid(shedBoundary, nodes, rectangles, boundaryNodes, resolution, internalPadding, raycastPadding))
	{
		LogMan::Log("Failed to generate FEM mesh", LOG_ERROR);
		ResetMeshs();
		return false;
	}
	
	ComputeBoundingBox(nodes, nodesSW, nodesNE);

	LogMan::Log("Succesfully generated mesh of " + std::to_string(rectangles.size()) + " rectangles!", LOG_SUCCESS);
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

	//look for triangles using this node and log warning if the change made it invalid
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		if (it->second.ContainsVertex(id))
		{
			if (!it->second.Validate())
				LogMan::Log("Warning! Updated triangle " + std::to_string(it->second.id) + " is invalid (Area: " + std::to_string(it->second.Area()) + ").", LOG_WARN);
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
	
	switch (params.meshType)
	{
	case ElementType::undefined:
		LogMan::Log("ERROR! Undefined mesh type. (state: Undefined)", LOG_ERROR);
		status = false;
		break;
	case ElementType::rectangle:
		if (rectangles.size() < 1) //existence of tris implies existence of nodes (though I probably need to protect both from change outside this file)
		{
			LogMan::Log("ERROR! No loaded mesh.", LOG_ERROR);
			status = false;
		}
		break;
	case ElementType::triangle:
		if (triangles.size() < 1) //existence of tris implies existence of nodes (though I probably need to protect both from change outside this file)
		{
			LogMan::Log("ERROR! No loaded mesh.", LOG_ERROR);
			status = false;
		}
		break;
	default:
		LogMan::Log("ERROR! Undefined mesh type (state: default).", LOG_ERROR);
		status = false;
		break;
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

void ConstructGlobalCapacitanceMatrix(ModelParameters const & params)
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
	
	//TODO replace this stupid switch with a unified loop once you figure out a way to abstract elements' containers.

	switch (params.meshType)
	{
		case ElementType::rectangle: //TODO consistant form is wrong for this
		{
			double diagVal, otherVal;
			for (auto it = rectangles.begin(); it != rectangles.end(); ++it)
			{
				//compute diagVal and otherVal
				if (params.useLumpedForm)
				{
					diagVal = it->second.Area() / 4.0; //(A/3.0) * 1.0
					otherVal = 0.0;
				}
				else
				{
					diagVal = 2.0 * it->second.Area() / 12.0;
					otherVal = it->second.Area() / 12.0; //(A/12.0) * 1.0
				}

				size_t vert[4] = { it->second.VertexID(0),it->second.VertexID(1),it->second.VertexID(2),it->second.VertexID(3) };

				globalC[vert[0]][vert[0]] += diagVal;
				globalC[vert[1]][vert[1]] += diagVal;
				globalC[vert[2]][vert[2]] += diagVal;
				globalC[vert[3]][vert[3]] += diagVal;

				globalC[vert[0]][vert[1]] += otherVal;
				globalC[vert[0]][vert[2]] += otherVal;
				globalC[vert[0]][vert[3]] += otherVal;

				globalC[vert[1]][vert[0]] += otherVal;
				globalC[vert[1]][vert[2]] += otherVal;
				globalC[vert[1]][vert[3]] += otherVal;

				globalC[vert[2]][vert[0]] += otherVal;
				globalC[vert[2]][vert[1]] += otherVal;
				globalC[vert[2]][vert[3]] += otherVal;

				globalC[vert[3]][vert[0]] += otherVal;
				globalC[vert[3]][vert[1]] += otherVal;
				globalC[vert[3]][vert[2]] += otherVal;
			}
		}
			break;

		case ElementType::triangle:
		{
			double diagVal, otherVal;
			for (auto it = triangles.begin(); it != triangles.end(); ++it)
			{
				//compute diagVal and otherVal
				if (params.useLumpedForm)
				{
					diagVal = it->second.Area() / 3.0; //(A/3.0) * 1.0
					otherVal = 0.0;
				}
				else
				{
					diagVal = 2.0 * it->second.Area() / 12.0;
					otherVal = it->second.Area() / 12.0; //(A/12.0) * 1.0
				}

				size_t const * vert = it->second.VertexIDs();

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
			break;

		default:
			break;
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

	Vector2D anchorNE(tiePoints[1][0] + (pixelScale[0] / 2.0), tiePoints[1][1] - (pixelScale[1] / 2.0));
	for (size_t row = 0; row < raster->Rows(); row++)
		for (size_t column = 0; column < raster->Columns(); column++)
		{
			Vector2D pixelPos(	anchorNE.x + column * pixelScale[0],
								anchorNE.y - row * pixelScale[1]);

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
double GetCurrentPrecipitation(double time, ModelParameters const & params, Element const & triangle)
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
		//size_t const * vert = it->second.vertIDs; //to simplify lines bellow.
		size_t vert[3] = { it->second.VertexID(0),it->second.VertexID(1), it->second.VertexID(2)};
		double newPrecipitation = GetCurrentPrecipitation(time + params.timeStep, params, it->second);
		double elementContrib = (it->second.Area() / 3.0) * (( 1.0 - params.femOmega) * it->second.elementPrecipitation + params.femOmega * newPrecipitation);
		elementContrib *= params.timeStep;
		it->second.elementPrecipitation = newPrecipitation;

		//test
		if (!it->second.ContainsVertex(17) && !it->second.ContainsVertex(18))
			elementContrib = 0.0;
		//end test

		result[vert[0]] += elementContrib;
		result[vert[1]] += elementContrib;
		result[vert[2]] += elementContrib;
	}

	return result;
}

std::pair<double, double> ComputeVelocityComponents(size_t nodeID, Vector_f64 waterElevation) //waterElevation = nodeElevation + head
{
	double head = waterElevation[nodeID] - nodeElevation[nodeID];
	
	double u = 3600.0 * sqrt(nodeSlopeX[nodeID]) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];
	double v = 3600.0 * sqrt(nodeSlopeY[nodeID]) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];

	//test
	/*int dir = lround(nodeFDR[nodeID]);
	double signX = (dir == 2 || dir == 3 || dir == 4) ? -1.0 : 1.0;
	double signY = (dir == 1 || dir == 2 || dir == 8) ? 1.0 : -1.0;
	u *= signX;
	v *= signY;*/
	//end test

	return std::pair<double, double>(u, v);
}

Matrix_f64 ComputeGlobalCoefficientsMatrix(ModelParameters const & params, Vector_f64 const & newHeads)
{
	//[C] + w*dt*[K]
	Matrix_f64 secondTerm(nodes.size(), nodes.size());

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		//						|	ui (yj - yk) + vi (xk - xj)		ui (yk - yi) + vi (xi - xk)		ui (yi - yj) + vi (xj - xi)	|
		//[K_element] = 1/6	*	|	uj (yj - yk) + vj (xk - xj)		uj (yk - yi) + vj (xi - xk)		uj (yi - yj) + vj (xj - xi)	|
		//						|	uk (yj - yk) + vk (xk - xj)		uk (yk - yi) + vk (xi - xk)		uk (yi - yj) + vk (xj - xi)	|

		Matrix_f64 kMat(3, 3);

		for (auto it = triangles.begin(); it != triangles.end(); ++it)
		{
			Vector2D const & i = it->second.Node(0);
			Vector2D const & j = it->second.Node(1);
			Vector2D const & k = it->second.Node(2);

			auto uvi = ComputeVelocityComponents(it->second.VertexID(0), newHeads);
			auto uvj = ComputeVelocityComponents(it->second.VertexID(1), newHeads);
			auto uvk = ComputeVelocityComponents(it->second.VertexID(2), newHeads);

			kMat[0][0] = uvi.first * (j.y - k.y) + uvi.second * (k.x - j.x);
			kMat[0][1] = uvi.first * (k.y - i.y) + uvi.second * (i.x - k.x);
			kMat[0][2] = uvi.first * (i.y - j.y) + uvi.second * (j.x - i.x);

			kMat[0][0] = uvj.first * (j.y - k.y) + uvj.second * (k.x - j.x);
			kMat[0][1] = uvj.first * (k.y - i.y) + uvj.second * (i.x - k.x);
			kMat[0][2] = uvj.first * (i.y - j.y) + uvj.second * (j.x - i.x);

			kMat[0][0] = uvk.first * (j.y - k.y) + uvk.second * (k.x - j.x);
			kMat[0][1] = uvk.first * (k.y - i.y) + uvk.second * (i.x - k.x);
			kMat[0][2] = uvk.first * (i.y - j.y) + uvk.second * (j.x - i.x);

			kMat *= params.timeStep * params.femOmega / 6.0;

			for (int row = 0; row < 3; row++)
				for (int column = 0; column < 3; column++)
					secondTerm[it->second.VertexID(row)][it->second.VertexID(column)] += kMat[row][column];
		}
	}

	return globalC + secondTerm;
}

Matrix_f64 ComputeGlobalConductanceMatrix(ModelParameters const & params)
{
	//[C] - dt*(1-w)*[K]

	Matrix_f64 secondTerm(nodes.size(), nodes.size());

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		//						|	ui (yj - yk) + vi (xk - xj)		ui (yk - yi) + vi (xi - xk)		ui (yi - yj) + vi (xj - xi)	|
		//[K_element] = 1/6	*	|	uj (yj - yk) + vj (xk - xj)		uj (yk - yi) + vj (xi - xk)		uj (yi - yj) + vj (xj - xi)	|
		//						|	uk (yj - yk) + vk (xk - xj)		uk (yk - yi) + vk (xi - xk)		uk (yi - yj) + vk (xj - xi)	|

		Matrix_f64 kMat(3, 3);

		for (auto it = triangles.begin(); it != triangles.end(); ++it)
		{
			Vector2D const & i = it->second.Node(0);
			Vector2D const & j = it->second.Node(1);
			Vector2D const & k = it->second.Node(2);

			auto uvi = ComputeVelocityComponents(it->second.VertexID(0), heads);
			auto uvj = ComputeVelocityComponents(it->second.VertexID(1), heads);
			auto uvk = ComputeVelocityComponents(it->second.VertexID(2), heads);

			kMat[0][0] = uvi.first * (j.y - k.y) + uvi.second * (k.x - j.x);
			kMat[0][1] = uvi.first * (k.y - i.y) + uvi.second * (i.x - k.x);
			kMat[0][2] = uvi.first * (i.y - j.y) + uvi.second * (j.x - i.x);

			kMat[0][0] = uvj.first * (j.y - k.y) + uvj.second * (k.x - j.x);
			kMat[0][1] = uvj.first * (k.y - i.y) + uvj.second * (i.x - k.x);
			kMat[0][2] = uvj.first * (i.y - j.y) + uvj.second * (j.x - i.x);

			kMat[0][0] = uvk.first * (j.y - k.y) + uvk.second * (k.x - j.x);
			kMat[0][1] = uvk.first * (k.y - i.y) + uvk.second * (i.x - k.x);
			kMat[0][2] = uvk.first * (i.y - j.y) + uvk.second * (j.x - i.x);

			kMat *= params.timeStep * (1.0 - params.femOmega) / 6.0;

			for (int row = 0; row < 3; row++)
				for (int column = 0; column < 3; column++)
					secondTerm[it->second.VertexID(row)][it->second.VertexID(column)] += kMat[row][column];
		}
	}

	return globalC - secondTerm;
}

#pragma region Test

bool IsBoundaryNode(int id)
{
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == id)
			return true;
	return false;
}

Vector_f64 _new_h, _precipComp, _RHS;

void TestShowValues(ModelParameters const & params)
{
	std::cout << "\n  id     |  head | newH |||  |||   u0   |    v0   |   u1    |   v1    |  precip  |  RHS\n";
	for (size_t i = 0; i < nodes.size(); i++)
	{
		auto uv0 = ComputeVelocityComponents(i, heads);
		auto uv1 = ComputeVelocityComponents(i, _new_h);

		std::cout << std::setw(4) <<  i << " " << (IsBoundaryNode(i) ? "[B]" : "   ") <<  " | " <<
			std::fixed << std::setprecision(2) << std::setw(5) << heads[i] << " | " << std::setw(5) << _new_h[i] <<
			"|||" << (heads[i] < _new_h[i] ? "UP" : (heads[i] > _new_h[i] ? "DN" : "--")) << "|||" <<
			std::fixed << std::setprecision(2) << std::setw(7) << uv0.first << " | " << std::setw(7) << uv0.second << " | " << 
			std::fixed << std::setprecision(2) << std::setw(7) << uv1.first << " | " << std::setw(7) << uv1.second << " | " <<
			std::fixed << std::setprecision(2) << std::setw(7) << _precipComp[i] << " | " <<
			std::fixed << std::setprecision(2) << std::setw(7) << _RHS[i] << std::endl;
	}

	std::cout << "\nConductance Matrix:\n";
	(globalC -  ComputeGlobalConductanceMatrix(params) ).DisplayOnCLI(1);

	std::cout << "\nCoefficients Matrix:\n";
	(ComputeGlobalCoefficientsMatrix(params, _new_h) - globalC).DisplayOnCLI(1);

	std::cout << std::endl;
}
#pragma endregion


void ComputeRHSVector(double time, ModelParameters const & params, Vector_f64 & outRHS)
{
	//[GlobalConductanceMat] * {h_0} + precipComponent

	outRHS = ComputeGlobalConductanceMatrix(params) * heads + ComputePreciptationVector(time, params);
	
	//test
	_RHS = outRHS;
	_precipComp = ComputePreciptationVector(time, params);

	TestShowValues(params);
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

	//Special consideration. Since the boundary node listing includes our exit node, we have to manually remove it.
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == params.outletNode)
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
	ConstructGlobalCapacitanceMatrix(params);

#pragma endregion

#pragma region Test
	std::cout << "\n===================================================\n";
	std::cout << "Global Capacitance";
	std::cout << "\n===================================================\n";
	globalC.DisplayOnCLI(0);
	
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

	nodeElevation = Vector_f64(nodes.size()); //test
	heads = nodeElevation;

	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop");
	
	while (time <= params.endTime)
	{
		std::cout << "\n\n===================================================\n";
		LogMan::Log("At T= " + std::to_string(time));
		std::cout << "===================================================\n";
	
		Vector_f64 newHeads = heads +Vector_f64(nodes.size(), 0.05);

		//internal loop
		for (size_t i = 0; i <= params.maxInternalIterations; i++)
		{
			_new_h = newHeads; //test

			Vector_f64 RHS;
			ComputeRHSVector(time, params, RHS);
			
			Vector_f64 residuals; //needed by solver
			Vector_f64 fixedNewH;

			Matrix_f64 coeffMat = ComputeGlobalCoefficientsMatrix(params, newHeads);

			//Adjust system for boundary cond
			//https://finite-element.github.io/7_boundary_conditions.html
			for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
			{
				RHS[*it] = nodeElevation[*it];

				for (size_t i = 0; i < coeffMat.Columns(); i++)
					coeffMat[*it][i] = 0.0;
				
				coeffMat[*it][*it] = 1.0;
			}

			if (!Solve(coeffMat, RHS, fixedNewH, residuals, params))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}

			//force computed newHeads to be positive (not sure about this)
			for (size_t j = 0; j < fixedNewH.Rows(); j++)
				fixedNewH[j] = Max(fixedNewH[j], nodeElevation[j]);
			
			std::cout << "current Internal Residual: " << std::fixed << std::setprecision(10) <<  (newHeads - fixedNewH).Magnitude() << std::endl;
			if ((newHeads - fixedNewH).Magnitude() <= params.internalResidualTreshold)
			{
				newHeads = fixedNewH;
				break;
			}
			else if (i >= params.maxInternalIterations)
				LogMan::Log("Reached maxInternalIterations without reaching appropriate h", LOG_WARN);
			
			newHeads = fixedNewH;
		}

		std::cout << "heads result:" << std::endl;
		for (size_t i = 0; i < nodes.size(); i++)
			std::cout << i << "\t:\t" << heads[i] << std::endl;


		heads = newHeads;
		time += params.timeStep;

		//test
		std::cout << "Enter to proceed to next step\n";
		std::cin.sync();
		std::cin.get();

	}

	return true;
}
