#include "SpatialDataModule.hpp"

Vector_f64 nodeSlope, nodeManning, nodeFDR;
Vector_f64 nodeElevation;

Matrix_f64 const * dem = NULL, *manningRaster = NULL, *slopes = NULL, *fdr = NULL;

int demID, manningRasterID, slopesID, fdrID;

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
		if (i < boundary.size() - 1)
			segment = std::pair<Vector2D const *, Vector2D const *>(&boundary[i], &boundary[i + 1]);
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

//points must be a copy. See usage of this function in GetNodeSamplingSubBoundary_Tri()
void ConvexHull(std::vector<Vector2D> points, std::vector<Vector2D> & outBoundaryPoints)
{
	outBoundaryPoints.clear();
	std::vector<Vector2D>::const_iterator curNodeIt = nodes.begin();

	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		if (it->x < curNodeIt->x)
			curNodeIt = it;
	}

	outBoundaryPoints.push_back(*curNodeIt);
	curNodeIt = nodes.begin();

	Vector2D currentPoint;

	while (true)
	{
		currentPoint = *curNodeIt;

		for (auto it = nodes.begin(); it != nodes.end(); ++it)
		{
			Vector2D vec1 = *it - outBoundaryPoints.back();
			Vector2D vec2 = currentPoint - *it;

			double crossProd = (vec1.y * vec2.x) - (vec1.x * vec2.y);

			if (crossProd < 0.0)
				currentPoint = *it;
		}

		if (currentPoint == outBoundaryPoints.front())
			break;

		outBoundaryPoints.push_back(currentPoint);

		++curNodeIt;
		if (curNodeIt == nodes.end())
			curNodeIt = nodes.begin();
	}
}

void GetNodeSamplingSubBoundary_Rect(size_t nodeID, std::vector<Vector2D> & outSubBoundary, bool sampleOutsideMesh) //sampleOutsideMesh is for boundary nodes.
{
	outSubBoundary.clear();
	for (auto it = rectangles.begin(); it != rectangles.end(); ++it)
	{
		if (it->second.ContainsVertex(nodeID))
			outSubBoundary.push_back(it->second.Centroid());
	}

	//The unordered_map impl doesn't guarantee order, but order is important here. Fortuntely, for rect elements (and for our current\
		impl) we know that they are an x-y alligned grid, and that every non-boundary node has exactly four subBoundary also alligned,\
		 so we can order them by computing bounds, then assigning on the correct order based on that.
	//Note that this is only for internal nodes. Boundary nodes (and outlet node) are handled in the else block.
	if (!IsBoundaryNode(nodeID) && outSubBoundary.size() > 3)
	{
		Vector2D boundSW, boundNE;
		ComputeBoundingBox(outSubBoundary, boundSW, boundNE);

		outSubBoundary.clear();
		outSubBoundary.push_back(boundSW);
		outSubBoundary.push_back(Vector2D(boundNE.x, boundSW.y));
		outSubBoundary.push_back(boundNE);
		outSubBoundary.push_back(Vector2D(boundSW.x, boundNE.y));
	}
	else //Assuming the meshing is working properly, all nodes are part of at least one element, so outSubBoundary should have min size of 1.
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

	if (!IsBoundaryNode(nodeID) && outSubBoundary.size() > 3) //Note: second test is for outlet node. This is not guaranteed for delauney tri.
	{
		//Due to ordering issues with unordered_maps ang mesh generators, the outSubBoundary above may not have its node ordered\
		in a way that creates perfect boundary of the subRegion we want. So we reorder to do so with a simple convex hull\
		generating algorithm.
		ConvexHull(outSubBoundary, outSubBoundary);
	}
	else
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
			Vector2D pixelPos(anchorNE.x + column * pixelScale[0],
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

	if (!SamplePixelsValues(nodeSamplingSubBound, raster, rasterID, values) || values.size() < 1)
	{
		LogMan::Log("Error in sampling pixel values for node: " + std::to_string(nodeID), LOG_ERROR);
		return false;
	}

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
		return std::pair<double, double>(DBL_MIN, DBL_MIN);
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

	if (slopeValues.size() < 1) //implies fdrValues.size() also < 1.
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

	//nodeSlopeX = Vector_f64(nodes.size());
	//nodeSlopeY = Vector_f64(nodes.size());
	nodeSlope = Vector_f32(nodes.size());
	nodeFDR = Vector_f64(nodes.size());

	for (size_t i = 0; i < nodes.size(); i++)
	{
		/*std::pair<double, double> newSlopes = SampleAggregatedSlopes(i, SpatialSamplingMethod::average, 0.000001);
		nodeSlopeX[i] = newSlopes.first;
		nodeSlopeY[i] = newSlopes.second;*/
		//std::cout << "sampling for node: " << i << std::endl;
		if (!SampleAggregatedPixels(i, SpatialSamplingMethod::average, slopes, slopesID, nodeSlope[i]) ||
			!SampleAggregatedPixels(i, SpatialSamplingMethod::majority, fdr, fdrID, nodeFDR[i]))
			LogMan::Log("WARNING! Failed to sample slope/FDR for node: " + std::to_string(i), LOG_WARN);

		//convert the slope from percentage to m/m
		nodeSlope[i] = nodeSlope[i] / 100.0;
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
