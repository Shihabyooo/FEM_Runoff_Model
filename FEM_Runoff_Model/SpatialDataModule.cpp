#include "SpatialDataModule.hpp"

Vector_f64 nodeSlope, nodeManning, nodeFDR;
Matrix_f64 const *manningRaster = NULL, *slopes = NULL, *fdr = NULL;

int manningRasterID, slopesID, fdrID;

bool LoadInputRasters(ModelParameters const & params)
{
	LogMan::Log("Loading rasters");
	//Note Variable Precipitation rasters are not loaded initially, but loaded during model execution. A pass should, however, be carried initially
	//to create the time column of the time series vs raster path (so whe can tell when to load which raster). Caveats of this approach is
	//that, due to inefficient GeoTIFF parser, this would signifcantly lower the performance.
	//An alternative would be to preload a timeseries for each node (loop over rasters for each pair of coords). Caveats is that this would
	//need memory equivalent to at least Matrix_f64[nodeCount][tsLength] plus a double[rasterCount], where tsLength = raster count if all raster
	//fall within range of simulation time.
	//A third alternative solution is to create these per-node timeseries, but cache them to disk instead. Performance would be slower than memory
	//load, but faster than on-demand raster load. Memory would be much less than memory load.

	bool status = true;

	status = status && FileIO::LoadRaster(params.slopesPath, &slopesID, &slopes);
	status = status && FileIO::LoadRaster(params.fdrPath, &fdrID, &fdr);

	if (params.variableManningCoefficients)
		status = status && FileIO::LoadRaster(params.manningCoefficientRasterPath, &manningRasterID, &manningRaster);

	if (!status) //error already logged with the function calls above.
	{
		FileIO::UnloadRaster(slopesID);
		if (params.variableManningCoefficients)
			FileIO::UnloadRaster(manningRasterID);
	}

	return status;
}

void UnloadAllRasters()
{
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
			continue;
		
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

//TODO this function breaks for some nodes. See synthetic watershed when sampling is not nearest neighbour
//points must be a copy. See usage of this function in GetNodeSamplingSubBoundary()
void ConvexHull(std::vector<Vector2D> points, std::vector<Vector2D> & outBoundaryPoints)
{
	outBoundaryPoints.clear();
	std::vector<Vector2D>::const_iterator curNodeIt = points.begin();

	for (auto it = points.begin(); it != points.end(); ++it)
	{
		if (it->x < curNodeIt->x)
			curNodeIt = it;
	}

	outBoundaryPoints.push_back(*curNodeIt);
	curNodeIt = points.begin();

	Vector2D currentPoint;

	while (true)
	{
		for (auto it = points.begin(); it != points.end(); ++it)
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
		if (curNodeIt == points.end())
			curNodeIt = points.begin();
	}
}

void GetNodeSamplingSubBoundary(size_t nodeID, std::vector<Vector2D> & outSubBoundary, bool sampleOutsideMesh) //sampleOutsideMesh is for boundary nodes.
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

		outSubBoundary.push_back(pos + Vector2D(avgDist, avgDist));
		outSubBoundary.push_back(pos + Vector2D(avgDist, -1.0 * avgDist));
		outSubBoundary.push_back(pos + Vector2D(-1.0 * avgDist, -1.0 * avgDist));
		outSubBoundary.push_back(pos + Vector2D(-1.0 * avgDist, avgDist));
	}
}

double FDR2Angle(int fdr) //return angle in radians. Assumes AgNPS fdr format
{
	//agnps fdr format: 1 = north, 2 = NE, 3 = East, 4 = SE, 5 = South, 6 = SW, 7 = West, 8 = NW
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
		return false; //Error already logged in GetRasterMappingParameters();
	
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
		return std::pair<Vector2Int, double>(Vector2Int(-1, -1), 0.0); //Error already logged in GetRasterMappingParameters();
	
	size_t _row = 0, _column = 0;
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
	if (method == SpatialSamplingMethod::nearest)
	{
		std::pair<Vector2Int, double> nearestPixel = SampleNearestPixel(nodes[nodeID], raster, rasterID);
		outValue = nearestPixel.second;
		return nearestPixel.first.x >= 0.0;
	}

	//construct a region of influence for node, bu connecting the centroid for all elements that this node is part of.
	std::vector<Vector2D> nodeSamplingSubBound;

	GetNodeSamplingSubBoundary(nodeID, nodeSamplingSubBound, true);

	std::vector<double> values;

	if (!SamplePixelsValues(nodeSamplingSubBound, raster, rasterID, values) || values.size() < 1)
	{
		LogMan::Log("Error in sampling pixel values for node: " + std::to_string(nodeID), LOG_ERROR);
		return false;
	}

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

bool CacheManningCoefficients(ModelParameters const & params)
{
	nodeManning = Vector_f64(nodes.size());

	size_t counter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		if (!params.variableManningCoefficients)
			nodeManning[counter] = params.fixedManningCoeffient;
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
	nodeSlope = Vector_f32(nodes.size());
	nodeFDR = Vector_f64(nodes.size());

	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (!SampleAggregatedPixels(i, SpatialSamplingMethod::nearest, slopes, slopesID, nodeSlope[i]) ||
			!SampleAggregatedPixels(i, SpatialSamplingMethod::nearest, fdr, fdrID, nodeFDR[i]))
			LogMan::Log("WARNING! Failed to sample slope/FDR for node: " + std::to_string(i), LOG_WARN);

		//convert the slope from percentage to m/m
		nodeSlope[i] = nodeSlope[i] / 100.0;
	}

	return true;
}