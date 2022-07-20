#pragma once
#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//TODO add a cleanup method to clear the allocated memory (superTriangles, rasters) when program closes.
ElementType activeMeshType = ElementType::undefined;
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
Vector_f64 nodeSlopeX, nodeSlopeY, nodeManning, nodeFDR; //TODO consider pre-computing sqrt(slope-x)/manning and sqrt(slope-y)/manning and caching them instead
//Vector_f64 nodeSlope, nodeManning, nodeFDR;

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

	/*LogMan::Log("Loaded data with bounds: "
				+ std::to_string(min.x) + ", " + std::to_string(min.y) + " and "
				+ std::to_string(max.x) + ", " + std::to_string(max.y));*/
}

void ResetMeshs()
{
	nodes.clear();
	triangles.clear();
	rectangles.clear();
	boundaryNodes.clear();
	
	activeMeshType = ElementType::undefined;
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

bool GenerateMesh(MeshGeneratorParameters & params)
{
	LogMan::Log("Attempting to generate FEM mesh");
	ResetMeshs();

	if (!MeshGen::ValidateParameters(params))
	{
		LogMan::Log("ERROR! Invalid meshing parameters!", LOG_ERROR);
		return false;
	}

	params.boundary = &shedBoundary;

	void * elementsContainerPtr = NULL;
	
	switch (params.meshType)
	{
	case ElementType::rectangle:
		elementsContainerPtr = static_cast<void *>(&rectangles);
		break;
	case ElementType::triangle:
		if (params.useCustomNodes)
		{
			if (!FileIO::LoadCoordinatePairsCSV(params.inNodesListPath, nodes))
				return false;
			params.inNodesList = &nodes;
		}
		elementsContainerPtr = static_cast<void *>(&triangles);
		break;
	}

	bool status = MeshGen::GenerateMesh(params, elementsContainerPtr, &nodes, &boundaryNodes);

	if (status)
	{
		ComputeBoundingBox(nodes, nodesSW, nodesNE);
		activeMeshType = params.meshType;
	}
	else
		ResetMeshs();

	return status;
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

bool IsBoundaryNode(size_t nodeID)
{
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == nodeID)
			return true;
	return false;
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
	
	switch (activeMeshType)
	{
	case ElementType::undefined:
		LogMan::Log("ERROR! No loaded mesh. (state: Undefined)", LOG_ERROR);
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

//Works for both mesh types
void ConstructLumpedCapacitanceMatrix()
{
	//Lumped Capacitance matrix for each (triangular) element is 3x3 matrix
	//[C_e] = A/3 * |	1	0	0	|
	//				|	0	1	0	|
	//				|	0	0	1	|
	//Where A is the area of element.
	//Same for diagnoal matrix for rectangular elements, but with a 4 x 4 matrix.

	size_t elementCount = 0;
	int nodesPerElement = 0;

	switch (activeMeshType)
	{
	case ElementType::rectangle:
		elementCount = rectangles.size();
		nodesPerElement = 4;
		break;
	case ElementType::triangle:
		elementCount = triangles.size();
		nodesPerElement = 3;
		break;
	}

	for (size_t i = 0; i < elementCount; i++)
	{
		Element * element = NULL;

		switch (activeMeshType)
		{
		case ElementType::rectangle:
			element = new Element(rectangles[i]);
			break;
		case ElementType::triangle:
			element = new Element(triangles[i]);
			break;
		}

		for (int intNodeID = 0; intNodeID < element->NodeCount(); intNodeID++)
			globalC[element->VertexID(intNodeID)][element->VertexID(intNodeID)] += element->Area() / 3.0;


		delete element;
	}
}

//TODO Implement
void ConstructConsistantCapacitanceMatrix_Rect()
{
	LogMan::Log("WARNING! Consistant Capacitance matrix for rectanguler elements not yet implemented.", LOG_WARN);
}

void ConstructConsistantCapacitanceMatrix_Tri()
{
	//Consistent Capacitance matrix for each element is 3x3 matrix
	//[C_e] = A/12 *	|	2	1	1	|
	//					|	1	2	1	|
	//					|	1	1	2	|
	//Where A is the area of element.
	
	double diagVal, otherVal;
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		
		diagVal = 2.0 * it->second.Area() / 12.0;
		otherVal = it->second.Area() / 12.0; //(A/12.0) * 1.0

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

void ConstructGlobalCapacitanceMatrix(ModelParameters const & params)
{
	globalC = Matrix_f64(nodes.size(), nodes.size());

	if (params.useLumpedForm)
		ConstructLumpedCapacitanceMatrix();
	else
	{
		switch (activeMeshType)
		{
		case ElementType::rectangle:
			ConstructConsistantCapacitanceMatrix_Rect();
			break;
		case ElementType::triangle:
			ConstructConsistantCapacitanceMatrix_Tri();
			break;
		}
	}
}

//TODO this function is a slightly modified version of similarily named function in GridMesh generator. Merge them.
bool IsPointInsideBoundary(Vector2D const & point, Vector2D const & raySource, std::vector<Vector2D> const & boundary)
{
	//https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
	//https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
	
	if (boundary.size() < 2)
		return false;
	
	size_t counter = 0;
	double dX = point.x - raySource.x;
	double dY = point.y - raySource.y;

	for (size_t i = 0; i < boundary.size(); i++)
	{
		std::pair<Vector2D const *, Vector2D const *> segment;
		if(i < boundary.size() - 1)
			segment = std::pair<Vector2D const *, Vector2D const *>(&boundary[i], &boundary [i+1]);
		else
			segment = std::pair<Vector2D const *, Vector2D const *>(&boundary[i], &boundary[0]);

		double dX2 = segment.first->x - segment.second->x;
		double dY2 = segment.first->y - segment.second->y;
		double denominator = dX * dY2 - dY * dX2;

		if (denominator == 0.0)
		{
			//std::cout << "!!!!!!!! Caught zero denominator!\n"; //test
			continue;
		}

		double dX3 = point.x - segment.first->x;
		double dY3 = point.y - segment.first->y;

		double t = dX3 * dY2 - dY3 * dX2;
		t = t / denominator;
		double u = dX3 * dY - dY3 * dX;
		u = u / denominator;

		if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0)
			counter++;
	}

	return counter % 2 != 0;
}

void GetNodeSamplingSubBoundary_Rect(size_t nodeID, std::vector<Vector2D> & outSubBoundary, bool sampleOutsideMesh) //sampleOutsideMesh is for boundary nodes.
{
	outSubBoundary.clear();
	for (auto it = rectangles.begin(); it != rectangles.end(); ++it)
	{
		if (it->second.ContainsVertex(nodeID))
			outSubBoundary.push_back(it->second.Centroid());
	}

	if (outSubBoundary.size() < 3) //Assuming the meshing is working properly, all nodes are part of at least one element, so outSubBoundary should have min size of 1.
	{
		Vector2D const & pos = nodes[nodeID];
		
		if (sampleOutsideMesh)
		{
			Vector2D delta = outSubBoundary[0] - pos;

			if (outSubBoundary.size() == 2)
				outSubBoundary.pop_back(); //because the computations bellow will add it again, but we don't want to complexify the code to check which one does.

			outSubBoundary.push_back(Vector2D(pos.x - delta.x, outSubBoundary[0].y));
			outSubBoundary.push_back(Vector2D(pos.x - delta.x, pos.y - delta.y));
			outSubBoundary.push_back(Vector2D(outSubBoundary[0].x, pos.y - delta.y));
		}
		else
		{
			if (outSubBoundary.size() == 2)
			{
				outSubBoundary.push_back(outSubBoundary[0].x == outSubBoundary[1].x ? Vector2D(nodes[nodeID].x, outSubBoundary[1].y) : Vector2D(outSubBoundary[1].x, nodes[nodeID].y));
				outSubBoundary.push_back(outSubBoundary[0].x == outSubBoundary[1].x ? Vector2D(nodes[nodeID].x, outSubBoundary[0].y) : Vector2D(outSubBoundary[0].x, nodes[nodeID].y));
			}
			else
			{
				outSubBoundary.push_back(Vector2D(outSubBoundary[0].x, pos.y));
				outSubBoundary.push_back(Vector2D(pos.x, pos.y));
				outSubBoundary.push_back(Vector2D(pos.x, outSubBoundary[0].y));
			}
		}
	}
}

void GetNodeSamplingSubBoundary_Tri(size_t nodeID, std::vector<Vector2D> & outSubBoundary, bool sampleOutsideMesh) //sampleOutsideMesh is for boundary nodes.
{
	outSubBoundary.clear();
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		if (it->second.ContainsVertex(nodeID))
			outSubBoundary.push_back(it->second.Centroid());
	}

	//TODO testing outSubBoundary.size() to determine whether node is boundary or not does not work in the case of Delauney generated elements. Handle this.
	if (outSubBoundary.size() < 3) 
	{
		Vector2D const & pos = nodes[nodeID];
		
		//TODO replace
		//This is a placeholder implementation
		//get average distance to current sample boundnodes.

		double avgDist = 0.0;
		for (auto it = outSubBoundary.begin(); it != outSubBoundary.end(); ++it)
			avgDist += nodes[nodeID].DistanceTo(*it);
		avgDist = avgDist / static_cast<double>(outSubBoundary.size());
		outSubBoundary.clear();

		/*if (sampleOutsideMesh)
		{*/
			outSubBoundary.push_back(pos + Vector2D(avgDist, avgDist));
			outSubBoundary.push_back(pos + Vector2D(avgDist, -1.0 * avgDist));
			outSubBoundary.push_back(pos + Vector2D(-1.0 * avgDist, -1.0 * avgDist));
			outSubBoundary.push_back(pos + Vector2D(-1.0 * avgDist, avgDist));
		/*}
		else
		{

		}*/
	}
}

double FDR2Angle(int fdr) //return angle in radians. Assumes AgNPS fdr format
{
	switch (fdr)
	{
	case 1: //90 degrees
		return 1.5707963267949;
	case 2: //45 deg
		return 0.785398163397448;
	case 3: //0 deg
		return 0.0;
	case 4: //315 deg
		return 5.49778714378214;
	case 5: //270 deg
		return 4.71238898038469;
	case 6: //225 deg
		return 3.92699081698724;
	case 7: //180 de
		return 3.14159265358979;
	case 8: //135 deg
		return 2.35619449019234;
	default: //shouldn't happen
		return 0.0;
	}
}

bool SamplePixelsValues(std::vector<Vector2D> const & samplingRegion, Matrix_f64 const * raster, int rasterID, std::vector<double> & outValues)
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
		return false;
	}

	Vector2D anchorNE(tiePoints[1][0] + (pixelScale[0] / 2.0), tiePoints[1][1] - (pixelScale[1] / 2.0));

	Vector2D subRegionSW, subRegionNE;
	ComputeBoundingBox(samplingRegion, subRegionSW, subRegionNE);
	Rect subRegionRect(subRegionSW, subRegionNE);

	for (size_t row = 0; row < raster->Rows(); row++)
		for (size_t column = 0; column < raster->Columns(); column++)
		{
			Vector2D pixelPos(anchorNE.x + column * pixelScale[0],
				anchorNE.y - row * pixelScale[1]);

			//To avoid having to ray cast for all pixels of raster, we do a simple bounds check by testing whether\
			pixel centroid falls within boundary, if true, then we proceed to the finer (and more expensive) raycast test.

			if (subRegionRect.ContainsInclusive(pixelPos))
			{
				Vector2D rayStart(subRegionSW.x - 10.0f, pixelPos.y);
				if (IsPointInsideBoundary(pixelPos, rayStart, samplingRegion))
					outValues.push_back(raster->GetValue(row, column));
			}
		}

	//cleanup memory
	delete[] pixelScale;
	delete[] tiePoints[0];
	delete[] tiePoints[1];
	delete[] tiePoints;

	return true;
}

bool ComputeSlopesComponents(size_t nodeID,
							std::vector<double> const & slopeValues,
							std::vector<double> const & fdrValues,
							std::vector<double> & outSlopesX,
							std::vector<double> & outSlopesY)
{

	outSlopesX.clear();
	outSlopesY.clear();

	for (size_t i = 0; i < slopeValues.size(); i++)
	{
		if (isnan(fdrValues[i]) || isnan(slopeValues[i]))
		{
			LogMan::Log("WARNING! Caught NaN values for FDR or Slopes for node: " + std::to_string(nodeID), LOG_WARN);
			continue;
		}
		
		int dir = lround(fdrValues[i]);

		if (dir < 1 || dir > 7)
		{
			LogMan::Log("WARNING! Caught invalid FDR code sampling for nodeID: " + std::to_string(nodeID), LOG_WARN);
			continue;
		}

		double angle = FDR2Angle(dir);

		outSlopesX.push_back(slopeValues[i] * cos(angle) / 100.0);
		outSlopesY.push_back(slopeValues[i] * sin(angle) / 100.0);
	}

	return true;
}

//TODO add NearestNeighbour as option for SpatialSamplingMethod and merge this with SampleAggregatedSlopes
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

bool SampleAggregatedPixels(size_t nodeID, SpatialSamplingMethod method, Matrix_f64 const * raster, int rasterID, double & outValue, double tolerance = -1.0) //tolereance required only for majority aggregation
{
	//construct a region of influence for node, bu connecting the centroid for all elements that this node is part of.
	std::vector<Vector2D> nodeSamplingSubBound;

	switch (activeMeshType)
	{
	case ElementType::triangle:
	{
		GetNodeSamplingSubBoundary_Tri(nodeID, nodeSamplingSubBound, true);
	}
	break;
	case ElementType::rectangle:
	{
		GetNodeSamplingSubBoundary_Rect(nodeID, nodeSamplingSubBound, true);
	}
	break;
	case ElementType::undefined: //should never reach this point in the code, but still.
		LogMan::Log("ERROR! Internal error. (Undefined meshtype in SampleSpatialAveragedPixels().", LOG_ERROR);
		return false;
	default:
		LogMan::Log("ERROR! Internal error. (Default meshtype in SampleSpatialAveragedPixels().", LOG_ERROR);
		return false;
	}

	std::vector<double> values;

	if (!SamplePixelsValues(nodeSamplingSubBound, raster, rasterID, values))
		return false;

	std::pair<double, double> result;
	switch (method)
	{
	case SpatialSamplingMethod::average:
		outValue = Average(values);
		break;
	case SpatialSamplingMethod::median:
		outValue = Median(values);
		break;
	case SpatialSamplingMethod::majority:
		outValue = Majority(values);
		break;
	default:
		LogMan::Log("ERROR! Internal error (Unspecified spatial sampling method in SampleAggregatedPixels())", LOG_ERROR);
		break;
	}

	return true;
}

//returns DBL_MIN values if error
std::pair<double, double> SampleAggregatedSlopes(size_t nodeID, SpatialSamplingMethod method, double tolerance = -1.0)
{
	//construct a region of influence for node, bu connecting the centroid for all elements that this node is part of.
	std::vector<Vector2D> nodeSamplingSubBound;

	switch (activeMeshType)
	{
	case ElementType::triangle:
	{
		GetNodeSamplingSubBoundary_Tri(nodeID, nodeSamplingSubBound, true);
	}
		break;
	case ElementType::rectangle:
	{
		GetNodeSamplingSubBoundary_Rect(nodeID, nodeSamplingSubBound, true);
	}
		break;
	case ElementType::undefined: //should never reach this point in the code, but still.
		LogMan::Log("ERROR! Internal error. (Undefined meshtype in SampleSpatialAveragedPixels().", LOG_ERROR);
		return std::pair<double, double> (DBL_MIN, DBL_MIN);
	default:
		LogMan::Log("ERROR! Internal error. (Default meshtype in SampleSpatialAveragedPixels().", LOG_ERROR);
		return std::pair<double, double>(DBL_MIN, DBL_MIN);
	}

	std::vector<double> slopeValues, fdrValues;

	if (!SamplePixelsValues(nodeSamplingSubBound, slopes, slopesID, slopeValues) ||
		!SamplePixelsValues(nodeSamplingSubBound, fdr, fdrID, fdrValues))
		return std::pair<double, double>(DBL_MIN, DBL_MIN);
	
	if (slopeValues.size() != fdrValues.size())
	{
		LogMan::Log("ERROR! Mismatched FDR and Slope sampling for nodeID: " + std::to_string(nodeID) + ". Check input rasters.", LOG_ERROR);
		return std::pair<double, double>(DBL_MIN, DBL_MIN);
	}

	std::vector<double> slopesX, slopesY;

	if (slopeValues.size() < 1 ) //implies fdrValues.size() also < 1.
	{
		//LogMan::Log("No pixel within subregion for node : " + std::to_string(nodeID) + ". Falling back to nearest neighbour sampling.");
		slopeValues.push_back(SampleNearestPixel(nodes[nodeID], slopes, slopesID).second);
		fdrValues.push_back(SampleNearestPixel(nodes[nodeID], fdr, fdrID).second);
	}
	
	if (!ComputeSlopesComponents(nodeID, slopeValues, fdrValues, slopesX, slopesY))
		return std::pair<double, double>(DBL_MIN, DBL_MIN);

	std::pair<double, double> result;
	switch (method)
	{
	case SpatialSamplingMethod::average:
		result.first = Average(slopesX);
		result.second = Average(slopesY);
		break;
	case SpatialSamplingMethod::median:
		result.first = Median(slopesX);
		result.second = Median(slopesY);
		break;
	case SpatialSamplingMethod::majority:
		result.first = Majority(slopesX, tolerance);
		result.second = Majority(slopesY, tolerance);
		break;
	default:
		LogMan::Log("ERROR! Internal error (Unspecified spatial sampling method in SampleAggregatedPixels())", LOG_ERROR);
		break;
	}

	//TODO remove this and fdr vector after removing references to the latter in the rest of the code
	nodeFDR[nodeID] = lround(Majority(fdrValues));

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
	//nodeSlope = Vector_f32(nodes.size());
	nodeFDR = Vector_f64(nodes.size());

	for (size_t i = 0; i < nodes.size(); i++)
	{
		std::pair<double, double> newSlopes = SampleAggregatedSlopes(i, SpatialSamplingMethod::average, 0.000001);
		nodeSlopeX[i] = newSlopes.first;
		nodeSlopeY[i] = newSlopes.second;

		//if (!SampleAggregatedPixels(i, SpatialSamplingMethod::average, slopes, slopesID, nodeSlope[i]) ||
		//	!SampleAggregatedPixels(i, SpatialSamplingMethod::majority, fdr, fdrID, nodeFDR[i]))
		//	LogMan::Log("WARNING! Failed to sample slope/FDR for node: " + std::to_string(i), LOG_WARN);

		////convert the slope from percentage to m/m
		//nodeSlope[i] = nodeSlope[i] / 100.0;
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

//Precipitation returned as meters per hour.
//current impl doesn't need triangle, but later it would.
double GetCurrentPrecipitation(double time, ModelParameters const & params, Element const * element)
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

	size_t elementCount = 0;
	int nodesPerElement = 0;

	switch (activeMeshType)
	{
	case ElementType::rectangle:
		elementCount = rectangles.size();
		nodesPerElement = 4;
		break;
	case ElementType::triangle:
		elementCount = triangles.size();
		nodesPerElement = 3;
		break;
	}

	for (size_t i = 0; i < elementCount; i++)
	{
		Element * element = NULL;

		switch (activeMeshType)
		{
		case ElementType::rectangle:
			element = new Element(rectangles[i]);
			break;
		case ElementType::triangle:
			element = new Element(triangles[i]);
			break;
		}

		double newPrecipitation = GetCurrentPrecipitation(time + params.timeStep, params, element);
		double elementContrib = (element->Area() / 3.0) * ((1.0 - params.femOmega) * element->elementPrecipitation + params.femOmega * newPrecipitation);
		elementContrib *= params.timeStep;
		element->elementPrecipitation = newPrecipitation;

		//test
		/*if (!element->ContainsVertex(1) && !element->ContainsVertex(2) && !element->ContainsVertex(3))
			elementContrib = 0.0;*/
		//end test

		for (int intNodeID = 0; intNodeID < element->NodeCount(); intNodeID++)
			result[element->VertexID(intNodeID)] += elementContrib;
		
		delete element;
	}

	return result;
}

//TODO have to vector_f64 to store computed nodal velocities, and a function to compute the newheads nodal velocities every \
iteration. newHead velos are then moved to oldHead velos after internal iterations are done.
std::pair<double, double> ComputeVelocityComponents(size_t nodeID, Vector_f64 waterElevation) //waterElevation = nodeElevation + head
{
	double head = waterElevation[nodeID] - nodeElevation[nodeID];
	
	/*double velocity = 3600.0 * sqrt(nodeSlope[nodeID]) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];

	double angle = FDR2Angle(lround(nodeFDR[nodeID]));

	double u = velocity * cos(angle);
	double v = velocity * sin(angle);*/

	double u = Sign(nodeSlopeX[nodeID]) * 3600.0 * sqrt(abs(nodeSlopeX[nodeID])) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];
	double v = Sign(nodeSlopeY[nodeID]) *3600.0 * sqrt(abs(nodeSlopeY[nodeID])) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];

	return std::pair<double, double>(u, v);
}

Matrix_f64 ElementKMatrix_Tri(Triangle const & tri, std::pair<double, double> uv[3])
{
	//						|	ui (yj - yk) + vi (xk - xj)		ui (yk - yi) + vi (xi - xk)		ui (yi - yj) + vi (xj - xi)	|
	//[K_element] = 1/6	*	|	uj (yj - yk) + vj (xk - xj)		uj (yk - yi) + vj (xi - xk)		uj (yi - yj) + vj (xj - xi)	|
	//						|	uk (yj - yk) + vk (xk - xj)		uk (yk - yi) + vk (xi - xk)		uk (yi - yj) + vk (xj - xi)	|

	//The 1/6 scalar is mulitplied in the calling function.

	Vector2D const & i = tri.Node(0);
	Vector2D const & j = tri.Node(1);
	Vector2D const & k = tri.Node(2);

	Matrix_f64 kMat(3, 3);

	//uv[0].second = uv[1].second = uv[2].second = (uv[0].second + uv[1].second + uv[2].second) / 3.0;//test
	//uv[0].first = uv[1].first = uv[2].first = (uv[0].first + uv[2].first + uv[2].first) / 3.0;//test

	kMat[0][0] = uv[0].first * (j.y - k.y) + uv[0].second * (k.x - j.x);
	kMat[0][1] = uv[0].first * (k.y - i.y) + uv[0].second * (i.x - k.x);
	kMat[0][2] = uv[0].first * (i.y - j.y) + uv[0].second * (j.x - i.x);

	kMat[1][0] = uv[1].first * (j.y - k.y) + uv[1].second * (k.x - j.x);
	kMat[1][1] = uv[1].first * (k.y - i.y) + uv[1].second * (i.x - k.x);
	kMat[1][2] = uv[1].first * (i.y - j.y) + uv[1].second * (j.x - i.x);

	kMat[2][0] = uv[2].first * (j.y - k.y) + uv[2].second * (k.x - j.x);
	kMat[2][1] = uv[2].first * (k.y - i.y) + uv[2].second * (i.x - k.x);
	kMat[2][2] = uv[2].first * (i.y - j.y) + uv[2].second * (j.x - i.x);

	return kMat;
}

Matrix_f64 ComputeGlobalCoefficientsMatrix(ModelParameters const & params, Vector_f64 const & newHeads)
{
	//[C] + w*dt*[K]
	Matrix_f64 coefMat(nodes.size(), nodes.size());
	Matrix_f64 kMat;
	Matrix_f64 cMat(3, 3);

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		std::pair<double, double> uv[3] = { ComputeVelocityComponents(it->second.VertexID(0), newHeads),
											ComputeVelocityComponents(it->second.VertexID(1), newHeads),
											ComputeVelocityComponents(it->second.VertexID(2), newHeads) };

		kMat = ElementKMatrix_Tri(it->second, uv) * (params.timeStep * params.femOmega / 6.0);

		cMat[0][0] = cMat[1][1] = cMat[2][2] = (it->second.Area() / 3.0);

		for (int row = 0; row < 3; row++)
			for (int column = 0; column < 3; column++)
				coefMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] + kMat[row][column]);
	}

	return coefMat;
}

Matrix_f64 ComputeGlobalConductanceMatrix(ModelParameters const & params)
{
	//[C] - dt*(1-w)*[K]

	Matrix_f64 condMat(nodes.size(), nodes.size());
	Matrix_f64 kMat;
	Matrix_f64 cMat(3, 3);

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		std::pair<double, double> uv[3] = { ComputeVelocityComponents(it->second.VertexID(0), heads),
											ComputeVelocityComponents(it->second.VertexID(1), heads),
											ComputeVelocityComponents(it->second.VertexID(2), heads) };

		kMat = ElementKMatrix_Tri(it->second, uv) * (params.timeStep * (1.0 - params.femOmega) / 6.0);

		cMat[0][0] = cMat[1][1] = cMat[2][2] = (it->second.Area() / 3.0);

		for (int row = 0; row < 3; row++)
			for (int column = 0; column < 3; column++)
				condMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] - kMat[row][column]);
	}
	
	return condMat;
}

#pragma region Test

Vector_f64 _new_h, _precipComp, _RHS;
std::vector<std::pair<double, double>> qTS;
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


	size_t targetRowToShow = 47;
	std::cout << "\nConductance Matrix:\n";
	Matrix_f64 _globalConduc = ComputeGlobalConductanceMatrix(params);
	//(globalC -  ComputeGlobalConductanceMatrix(params) ).DisplayOnCLI(0);
	//ComputeGlobalConductanceMatrix(params).DisplayOnCLI(0);
	for (size_t i = 0; i < _globalConduc.Columns(); i++)
		std::cout << std::fixed << std::setprecision(0) << std::setw(2) << _globalConduc[targetRowToShow][i] << " ";
	std::cout << "\n";


	std::cout << "\nCoefficients Matrix:\n";
	Matrix_f64 _globalCoef = ComputeGlobalCoefficientsMatrix(params, _new_h);
	//(ComputeGlobalCoefficientsMatrix(params, _new_h) - globalC).DisplayOnCLI(0);
	//ComputeGlobalCoefficientsMatrix(params, _new_h).DisplayOnCLI(0);
	for (size_t i = 0; i < _globalCoef.Columns(); i++)
		std::cout << std::fixed << std::setprecision(0) << std::setw(2) << _globalCoef[targetRowToShow][i] << " ";
	std::cout << "\n";


	std::cout << std::endl;
}
#pragma endregion

void ComputeRHSVector(double time, ModelParameters const & params, Vector_f64 & outRHS)
{
	//[GlobalConductanceMat] * {h_0} + precipComponent

	outRHS = ComputeGlobalConductanceMatrix(params) * heads;// +ComputePreciptationVector(time, params);
	
	//test
	_RHS = outRHS;
	_precipComp = ComputePreciptationVector(time, params);

	//TestShowValues(params);
}

//[A]{x} = {b}
void AdjustForBoundaryConditions(Matrix_f64 & aMat, Vector_f64 & xVec, Vector_f64 & bVec)
{
	//Using https://finite-element.github.io/7_boundary_conditions.html
	/*for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	{
		bVec[*it] = nodeElevation[*it];

		for (size_t i = 0; i < aMat.Columns(); i++)
			aMat[*it][i] = 0.0;

		aMat[*it][*it] = 1.0;
	}*/

	//Using Istok's method
	size_t reducedSystemSize = nodes.size() - boundaryNodes.size();
	Matrix_f64 adjustedAMat(reducedSystemSize, reducedSystemSize);
	Vector_f64 adjustedXVec(reducedSystemSize);
	Vector_f64 adjustedBVec(reducedSystemSize);

	size_t rowCounter = 0;
	size_t columnCounter = 0;
	for (size_t row = 0; row < aMat.Rows(); row++)
	{
		if (IsBoundaryNode(row))
			continue;
		
		adjustedBVec[rowCounter] = bVec[row];
		//adjustedXVec[rowCounter] = xVec[row]; //pointless. xVec is not filled with data yet.

		for (size_t column = 0; column < aMat.Columns(); column++)
		{
			if (IsBoundaryNode(column))
				continue;

			adjustedAMat[rowCounter][columnCounter] = aMat[row][column];
			columnCounter++;
		}
		columnCounter = 0;
		rowCounter++;
	}
	
	aMat = adjustedAMat;
	bVec = adjustedBVec;
	xVec = adjustedXVec;
}

void ReIntroduceBoundaryNodes(Vector_f64 & vec)
{
	if (vec.Rows() == nodes.size()) 
		return;

	Vector_f64 adjustedVec(nodes.size());

	size_t rowCounter = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (IsBoundaryNode(i))
		{
			adjustedVec[i] = 0.0;
		}
		else
		{
			adjustedVec[i] = vec[rowCounter];
			rowCounter++;
		}
	}

	vec = adjustedVec;
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
	/*std::cout << "\n===================================================\n";
	std::cout << "Global Capacitance";
	std::cout << "\n===================================================\n";
	globalC.DisplayOnCLI(0);*/
	
	std::cout << "\n===================================================\n";
	std::cout << "Time series";
	std::cout << "\n===================================================\n";
	for (int i = 0; i < params.unitTimeSeries.size; i++)
		std::cout << params.unitTimeSeries.series[i].first << " - " << params.unitTimeSeries.series[i].second << std::endl;
	std::cout << "\n===================================================\n";
	std::cout << "Nodes";
	std::cout << "\n===================================================\n";
	size_t nodeIDCounter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		std::cout << nodeIDCounter << "\t - \t" << std::fixed;
		Print(*it);
		nodeIDCounter++;
	}

	std::cout << "\n===================================================\n";
	std::cout << "Boundery Nodes";
	std::cout << "\n===================================================\n";
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		std::cout << *it << std::endl;
	
	std::cout << "\n===================================================\n";
	std::cout << "Elevations, Slopes, Manning roughness coef";
	std::cout << "\n===================================================\n";


	std::cout << "node | Elev  |   n  |  Sx |  Sy |  FDR  \n";
	for (size_t i = 0; i < nodes.size(); i++)
		std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " << nodeElevation[i] << " | " <<
		nodeManning[i] << " | " << nodeSlopeX[i] << " | " << nodeSlopeY[i] << " | "  << nodeFDR[i] << std::endl;
	/*std::cout << "node | Elev  |   n  |  Slope  |  FDR  \n";
	for (size_t i = 0; i < nodes.size(); i++)
		std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " << nodeElevation[i] << " | "  <<
					nodeManning[i] << " | " << nodeSlope[i] << " | " << nodeFDR[i]<< std::endl;*/
	
	//return false;
#pragma endregion

	nodeElevation = Vector_f64(nodes.size()); //test
	heads = nodeElevation;

	heads[57] += 2.0;//test

	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop");
	
	while (time <= params.endTime)
	{
		std::cout << "\n\n===================================================\n";
		LogMan::Log("At T= " + std::to_string(time));
		std::cout << "===================================================\n";
	
		//Vector_f64 newHeads = heads + Vector_f64(nodes.size(), 0.001);
		Vector_f64 newHeads(nodes.size());
		for (size_t i = 0; i < nodes.size(); i++)
			if (!IsBoundaryNode(i))
				newHeads[i] = heads[i] + 0.001;


		//internal loop
		for (size_t i = 0; i <= params.maxInternalIterations; i++)
		{
			_new_h = newHeads; //test

			Vector_f64 RHS;
			ComputeRHSVector(time, params, RHS);
			
			Vector_f64 residuals; //needed by solver
			Vector_f64 fixedNewH;

			Matrix_f64 coeffMat = ComputeGlobalCoefficientsMatrix(params, newHeads);

			AdjustForBoundaryConditions(coeffMat, fixedNewH, RHS);

			/*if (!Solve(coeffMat, RHS, fixedNewH, residuals, params))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}*/

			//test
			if (!SolverSOR(coeffMat, RHS, fixedNewH, residuals, 0.75))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}

			//force computed newHeads to be positive (not sure about this)
			for (size_t j = 0; j < fixedNewH.Rows(); j++)
				fixedNewH[j] = Max(fixedNewH[j], nodeElevation[j]);
			
			ReIntroduceBoundaryNodes(fixedNewH);

			std::cout << "Solver residual: " << std::fixed << std::setprecision(10) << residuals.Magnitude() << std::endl; //test
			std::cout << "current Internal Residual: " << std::fixed << std::setprecision(10) <<  (newHeads - fixedNewH).Magnitude() << std::endl; //test
			if ((newHeads - fixedNewH).Magnitude() <= params.internalResidualTreshold)
			{
				newHeads = fixedNewH;
				break;
			}
			else if (i >= params.maxInternalIterations)
				LogMan::Log("Reached maxInternalIterations without reaching appropriate h", LOG_WARN);
			
			newHeads = fixedNewH;
		}
		TestShowValues(params);

		std::cout << "\n------------------------------------------------------\n";
		std::cout << "heads result at time: " << std::setprecision(4) << time << " --> " << std::setprecision(4) << time + params.timeStep << std::endl;
		for (size_t i = 0; i < heads.Rows(); i++)
			std::cout << i << "\t slopes: " << std::setprecision(4) << nodeSlopeX[i] << " - " << nodeSlopeY[i]  << "\t - \t" << std::setprecision(4) << heads[i] << "\t-->\t" <<
						std::setprecision(4) << newHeads[i] << "  " <<
						(newHeads[i] > heads[i] ? "UP" : (newHeads[i] == heads[i] ? "--" : "DN"))<<
						std::endl;
		std::cout << "\n------------------------------------------------------\n";
		double qx = sqrt(nodeSlopeX[params.outletNode]) * pow(heads[params.outletNode], 5.0 / 3.0) / nodeManning[params.outletNode];
		double qy = sqrt(nodeSlopeY[params.outletNode]) * pow(heads[params.outletNode], 5.0 / 3.0) / nodeManning[params.outletNode];
		double q = sqrt(qx * qx + qy * qy);
		double flowWidth = 2.0 *  abs((triangles[0].Centroid() - triangles[0].Node(0)).x); //test. 
		std::cout << "h: " << std::setprecision(7) << heads[params.outletNode]  << "\tQ: "  << std::setprecision(3) << (q * flowWidth) << std::endl;
		qTS.push_back(std::pair<double, double>(time, q * flowWidth));
		std::cout << "\n------------------------------------------------------\n";

		heads = newHeads;
		time += params.timeStep;

		//test
		std::cout << "Enter to proceed to next step\n";
		std::cin.sync();
		std::cin.get();
	}

	//test
	std::cout << " time \t Q (cms)\n";
	for (auto it = qTS.begin(); it != qTS.end(); ++it)
		std::cout << std::setw(6) << std::setfill('0') << std::setprecision(3) << it->first << "\t" << std::setprecision(5) << it->second << std::endl;
	qTS.clear();

	return true;
}
