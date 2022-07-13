#include "GridMesh.hpp"

//Clear the output vectors/unord_maps.
//Compute absolute grid values
	//vertical number of cells = resolution, rows = resolution + 1 (since the rows are for the nodes, not elements)
	//the element height = (delta.y - 2.0 * padding) / resolution
	//y pos of first row = boundSE.y + padding. then add element.Height for every following row.
	//take element width = height.
	//number of (element) columns =  floor( (delta.x - 2 * padding) / resolution).
	//x pos of first column = boundSE.x + padding, then add element.width for every following column
//bin the line segments of boundary as follows:
	//Create an array of vectors of size = rows.
	//loop over rows, store row y coordinate
		//loop over segments. e.g. for (size_t i = 0; i < boundary.size(); i++)\
			we create a segment as pair<vector2d*, vector2d*> = {&boundary[i], &boundary[i+1]}. With exception when we hit last point\
			we loopback to {&boundary[i], &boundary[0]}.
				//test that the row's y is between the segment's begin y and end y
				//e.g. if ((row.y - begin.y) * (row.y - end.y)) < 0.0) -> row is between segment ends.
					//if true, push segment pair appropriate row's vector in the array above.
//The binning step can be skipped for shorter code, but then we'll have to raycast test all elements (or at least do the\
check above) for every point tested bellow. With the binning, we can limit the internal loop tests to only binned segments.
//
//Using number of rows and columns above, creat a 2d array of bools, this will for cheaper mapping of points within the mesh\
and later generation of elements. All bools are initialized to false (points are not within mesh).
//Loop over each row
	//loop over each column (now we have a point to test)
		//Test that this point is inside boundary by raycasting
		//for the point, the ray (more like line) to is cast from outside the boundary to the left (startPosition.x = boundarySE.x - rayPadding.\
		y is the same. rayPadding should be large compared to padding), towards the point.
		//start a counter at zero.
		//loop over the binned segments for current row
			//test interesction between segment and ray
				//if interescts, counter++
		//check counter%2, if 0 point is outside boundary, if 1, point is inside.
			//for inside points, mark the corrosponding bool in the map as true.
//After the end of the loop we should have  a list of nodes that all are within the boundary, but not all are usefull for elements (e.g.\
a point may have no adjacent node and equivalent pair of nodes above to create an element from.
//create a map of size_ts of same dimensions, this will store node IDs
//start counter = 0;
//Loop again over rows and columns and test the bool map
	//if this point is true
		//test d8 adjacent points in groups of threes (four quadrants) to a single bool.\
		i.e bool q1 = N && NE && E; q2 = E && SE && S, q3 = S && SW && W, q4 = W && NW && N. \
		The value of this node would then be q1 || q2 || q3 || q\
		Note that for boundary cases (nodes at end of the map/grid), the direction equivalent to boundary is false.\
		This test would ensure that every point in the map can be used to make an element with adjacent nodes.
		//If the value of the node is true:
			//id = counter; counter++;
			//compute the coordinates of the node (using parameters above) and push to nodesList.
//At this point, we have a map of nodes that can be used to create elements.
//loop over rows and columns
	//check that point is true.
		//check that it is a SW corner (i.e. nodes to E, N and NE are true).
		//if true
			//create an element from this node, an N, E, NE nodes, push to elementsList.
//Free any temp memory created above.


//NOTE! IsNodeValid is also used for determining whether a node is boundary node or not. Simply put, a boundary node is one\
that has at least one of its quadrants without an element. An internal node is surrounded by elements on all quadrants.

struct GriddingParameters
{
public:

	bool Validate()
	{
		if (elementWidth < MIN_RECTANGLE_WIDTH || elementHeight < MIN_RECTANGLE_HEIGHT)
		{
			LogMan::Log("ERROR! Resulting element width/height are too small/invalid at: "
				+ std::to_string(elementWidth) + "x" + std::to_string(elementHeight), LOG_ERROR);
		
			return false;
		}
		
		if (nodesRows < 2 || nodesColumns < 2)
		{
			LogMan::Log("ERROR! Resulting nodes rows/columns are too few at: "
				+ std::to_string(nodesRows) + "x" + std::to_string(nodesColumns), LOG_ERROR);

			return false;
		}

		return true;
	}

	double elementWidth = MIN_RECTANGLE_WIDTH;
	double elementHeight = MIN_RECTANGLE_HEIGHT;
	size_t nodesRows = 2;
	size_t nodesColumns = 2;
	Vector2D gridAnchorSW = Vector2D();
};

std::vector<std::pair<Vector2D const *, Vector2D const *>> * binnedSegments = NULL;
bool ** nodesMap = NULL;
size_t ** idsMap = NULL;

GriddingParameters params;

inline Vector2D GridToWorldPos(size_t row, size_t column)
{
	return Vector2D(params.gridAnchorSW.x + column * params.elementWidth,
					params.gridAnchorSW.y + row * params.elementHeight);
}

bool ValidateInput(std::vector<Vector2D> const & boundary, size_t resolution, double internalPadding, double rayCastPadding)
{
	bool status = true;
	if (boundary.size() < MIN_BOUNDARY_TO_DISCRETIZE)
	{
		LogMan::Log("ERROR! Cannot work with boundry with smaller edges than " + std::to_string(MIN_BOUNDARY_TO_DISCRETIZE), LOG_ERROR);
		status = false;
	}
	
	if (resolution < MIN_RESOLUTION)
	{
		LogMan::Log("ERROR! Resolution must be >= " + std::to_string(MIN_RESOLUTION), LOG_ERROR);
		status = false;
	}
	
	if (internalPadding <= 0.0)
	{
		LogMan::Log("ERROR! Internal padding must be greater than 0.0" , LOG_ERROR);
		status = false;
	}
	
	if (rayCastPadding <= 0.0)
	{
		LogMan::Log("ERROR! Raycast padding must be greater than 0.0", LOG_ERROR);
		status = false;
	}

	return status;
}

void BoundingBox(std::vector<Vector2D> const & points, Vector2D & min, Vector2D & max)
{
	//TODO I really should centralize this...
	min = points[0];
	max = points[0];

	for (auto it = points.begin(); it < points.end(); it++)
	{
		min.x = Min(it->x, min.x);
		min.y = Min(it->y, min.y);
		max.x = Max(it->x, max.x);
		max.y = Max(it->y, max.y);
	}
}

bool ComputeGriddingParameters(std::vector<Vector2D> const & boundary, double internalPadding, double resolution)
{
	Vector2D sw, ne;
	BoundingBox(boundary, sw, ne);
	Vector2D delta = ne - sw;

	if ( delta.x < MIN_RECTANGLE_WIDTH || delta.y < MIN_RECTANGLE_HEIGHT) //TOD refine this check.
	{
		LogMan::Log("ERROR! Supplied bad boundary to GridMesh (Failed grid param compute stage)", LOG_ERROR);
		return false;
	}
	
	params.elementWidth = params.elementHeight = (delta.y - 2.0 * internalPadding) / static_cast<double>(resolution);
	params.gridAnchorSW = sw + Vector2D(internalPadding, internalPadding);
	params.nodesRows = resolution + 1;
	params.nodesColumns = llround(floor(delta.x - 2.0 * internalPadding) / static_cast<double>(params.elementWidth));

	//test
	std::cout << "gridding parameters\n";
	std::cout << "element size: " << params.elementWidth << " x " << params.elementHeight << std::endl;
	std::cout << "rows x columns: " << params.nodesRows << " x " << params.nodesColumns << std::endl;
	std::cout << "anchor: ";
	Print(params.gridAnchorSW);
	//endtest

	return params.Validate();
}

bool BinBoundarySegments(std::vector<Vector2D> const & boundary)
{
	if (binnedSegments != NULL)
		delete[] binnedSegments;

	binnedSegments = new std::vector<std::pair<Vector2D const *, Vector2D const *>>[params.nodesRows];

	for (size_t row = 0; row < params.nodesRows; row++)
	{
		Vector2D pointPos = GridToWorldPos(row, 0);
		for (size_t i = 0; i < boundary.size(); i++)
		{
			std::pair<Vector2D const *, Vector2D const *> segment;
				
			if (i < boundary.size() - 1) //first to second-to-last segements
				segment = std::pair<Vector2D const *, Vector2D const *>(&boundary[i], &boundary[i+1]);
			else //last segment
				segment = std::pair<Vector2D const *, Vector2D const *>(&boundary[i], &boundary[0]);
				
			double testValue = (pointPos.y - segment.first->y) * (pointPos.y - segment.second->y);
			if (testValue <= 0.0) //TODO check whether the intersection test works when the intersection is right at the end.
				binnedSegments[row].push_back(segment);
		}
		if (binnedSegments[row].size() < 2) //Though this shouldn't happen...
		{
			LogMan::Log("ERROR! Supplied bad boundary to GridMesh (Failed binning stage)", LOG_ERROR);
			return false;
		}
		//std::cout << "Binned at row: " << row << " = " << binnedSegments[row].size() << std::endl; //test
	}
	return true;
}

bool AllocateNodesMaps()
{
	try
	{
		nodesMap = new bool*[params.nodesRows];
		bool * helperPtr = new bool[params.nodesRows * params.nodesColumns]();

		idsMap = new size_t *[params.nodesRows];
		size_t * helperPtr2 = new size_t[params.nodesRows * params.nodesColumns]();

		for (size_t i = 0; i < params.nodesRows; i++)
		{
			nodesMap[i] = helperPtr;
			idsMap[i] = helperPtr2;

			helperPtr += params.nodesColumns;
			helperPtr2 += params.nodesColumns;
		}

	}
	catch (const std::bad_alloc& exception)
	{
		LogMan::Log("ERROR! Could not allocate internal map. " + std::string(exception.what()), LOG_ERROR);
		return false;
	}

	return true;
}

bool IsPointInsideBoundary(Vector2D const & point, Vector2D const & raySource, size_t row)
{
	//https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
	//https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
	size_t counter = 0;
	double dX = point.x - raySource.x;
	double dY = point.y - raySource.y;
	
	//std::cout << "testing point: ";//test
	//Print(point);//test

	for (size_t i = 0; i < binnedSegments[row].size(); i++)
	{
		std::pair<Vector2D const *, Vector2D const *> const & segment = binnedSegments[row][i];

		//std::cout << "testing segment: "; //test
		//Print(*segment.first, true); //test
		//std::cout << " | "; //test
		//Print(*segment.second); //test


		double dX2 = segment.first->x - segment.second->x;
		double dY2 = segment.first->y - segment.second->y;
		double denominator = dX * dY2 - dY * dX2;
		
		if (denominator == 0.0)
		{
			std::cout << "!!!!!!!! Caught zero denominator!\n"; //test
			continue;
		}

		double dX3 = point.x - segment.first->x;
		double dY3 = point.y - segment.first->y;

		double t = dX3 * dY2 - dY3 * dX2;
		t = t / denominator;
		double u = dX3 * dY - dY3 * dX;
		u = u / denominator;

		//std::cout << "T: " << t << ", u: " << u << std::endl;

		if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0)
			counter++;
	}

	return counter % 2 != 0;
}

//Validates node, if valid, assign an id to it, if is boundary, push id to outBoundaryNodes.
bool ValidateNode(size_t row, size_t column, size_t idToAssign, std::vector<size_t> & outBoundaryNodes) //checks that node  creates a valid rect in at least one of its four quadrants
{
	//test d8 adjacent points in groups of threes (four quadrants) to a single bool.\
		i.e bool q1 = N && NE && E; q2 = E && SE && S, q3 = S && SW && W, q4 = W && NW && N. \
		The value of this node would then be q1 || q2 || q3 || q\
		Note that for boundary cases (nodes at end of the map/grid), the direction equivalent to boundary is false.\
		This test would ensure that every point in the map can be used to make an element with adjacent nodes.
		//If the value of the node is true:
	if (!nodesMap[row][column]) //Even though this should be enforced in the calling code...
		return false;

	bool q1 = false, q2 = false, q3 = false, q4 = false;

	if (row < params.nodesRows - 1)
	{
		if (column < params.nodesColumns - 1)
			q1 = nodesMap[row][column + 1] && nodesMap[row + 1][column] && nodesMap[row + 1][column + 1]; 
		if (column > 0)
			q4 = nodesMap[row][column - 1] && nodesMap[row + 1][column] && nodesMap[row + 1][column - 1];
	}
	if (row > 0)
	{
		if (column < params.nodesColumns - 1)
			q2 = nodesMap[row][column + 1] && nodesMap[row - 1][column] && nodesMap[row - 1][column + 1];
		if (column > 0)
			q3 = nodesMap[row][column - 1] && nodesMap[row - 1][column] && nodesMap[row - 1][column - 1];
	}

	bool isValid = q1 || q2 || q3 || q4;
	bool isInternal = q1 && q2 && q3 && q4;

	if (isValid)
	{
		idsMap[row][column] = idToAssign;
		if (!isInternal)
			outBoundaryNodes.push_back(idToAssign);
	}
	
	//Warning for potentially problematic nodes (not d4 connected, i.e. doesn't create shared edges)
	if ((!(q1 || q3) && (q2 && q4)) || (!(q2 || q4) && (q1 && q3)))
		LogMan::Log("Warning! D4-discontinuity at node " + std::to_string(idToAssign), LOG_WARN);
		
	return isValid;
}

bool IsSWNode(size_t row, size_t column, size_t ** outRectNodes) //outRectNodes is size_t[4]
{
	//check not at eastern nor northern edge
	if (row >= params.nodesRows - 1 || column >= params.nodesColumns - 1)
		return false;
	
	(*outRectNodes)[0] = idsMap[row][column];
	(*outRectNodes)[1] = idsMap[row][column + 1];
	(*outRectNodes)[2] = idsMap[row + 1][column];
	(*outRectNodes)[3] = idsMap[row + 1][column + 1];

	//check nodes at [row][column + 1], [row + 1][column] and [row + 1][column + 1] are all true.
	return nodesMap[row][column + 1] && nodesMap[row + 1][column] && nodesMap[row + 1][column + 1];
}

bool GenerateNodes(std::vector<Vector2D> const & boundary, std::vector<Vector2D> & outNodes, std::vector<size_t> & outBoundaryNodes, double rayCastPadding)
{
	LogMan::Log("Generating nodes.");

	if (!AllocateNodesMaps())
		return false;

	//first pass
	for (size_t row = 0; row < params.nodesRows; row++)
	{
		for (size_t column = 0; column < params.nodesColumns; column++)
		{
			Vector2D pointPos = GridToWorldPos(row, column);
			Vector2D raySource(params.gridAnchorSW.x - rayCastPadding, pointPos.y);
			if (IsPointInsideBoundary(pointPos, raySource, row))
				nodesMap[row][column] = true;
		}
	}

	//second pass
	size_t counter = 0;
	for (size_t row = 0; row < params.nodesRows; row++)
	{
		for (size_t column = 0; column < params.nodesColumns; column++)
		{
			if (nodesMap[row][column])
			{
				if (ValidateNode(row, column, counter, outBoundaryNodes))
				{
					//idsMap[row][column] = counter;
					outNodes.push_back(GridToWorldPos(row, column));
					counter++;
				}
				else
					nodesMap[row][column] = false;
			}
		}
	}

	LogMan::Log("Generated " + std::to_string(outNodes.size()) + " nodes.");
	LogMan::Log("With " + std::to_string(outBoundaryNodes.size()) + " boundary nodes.");
	return true;
}

void GenerateElements(std::vector<Vector2D> const & nodes, std::unordered_map<size_t, Rectangle> & outRectList)
{
	//loop over rows and columns
	//check that point is true.
		//check that it is a SW corner (i.e. nodes to E, N and NE are true).
		//if true
			//create an element from this node, an N, E, NE nodes, push to elementsList.

	LogMan::Log("Attempting to generate elements");

	size_t counter = 0;
	size_t * rectNodes = new size_t[4];

	for (size_t row = 0; row < params.nodesRows; row++)
	{
		for (size_t column = 0; column < params.nodesColumns; column++)
		{
			if (nodesMap[row][column] && IsSWNode(row, column, &rectNodes))
			{
				//std::cout << "creating rect of nodes: " << rectNodes[0] << ", " << rectNodes[1] << ", " << rectNodes[2] << ", " << rectNodes[3] << std::endl; //test
				Rectangle newRect(counter, rectNodes, &nodes);
				outRectList.insert({ counter, newRect });
				counter++;
			}
		}
	}

	delete[] rectNodes;

	LogMan::Log("Generated " + std::to_string(outRectList.size()) + " elements.");
}

void MemoryCleanup()
{
	if (binnedSegments != NULL)
		delete[] binnedSegments;
	binnedSegments = NULL;

	if (nodesMap != NULL)
	{
		delete[] nodesMap[0];
		delete[] nodesMap;
		nodesMap = NULL;
	}

	if (idsMap != NULL)
	{
		delete[] idsMap[0];
		delete[] idsMap;
		idsMap = NULL;
	}
}

bool GenerateGrid(std::vector<Vector2D> const & boundary,
				std::vector<Vector2D> & outNodes,
				std::unordered_map<size_t, Rectangle> & outRectList,
				std::vector<size_t> & outBoundaryNodes,
				size_t resolution,
				double internalPadding, double rayCastPadding)
{
	//checks
	if (!ValidateInput(boundary, resolution, internalPadding, rayCastPadding))
		return false;

	LogMan::Log("Starting grid discretization.");

	outNodes.clear();
	outRectList.clear();
	outBoundaryNodes.clear();

	if (!ComputeGriddingParameters(boundary, internalPadding, resolution))
		return false;

	//from now on, we have to call MemoryCleanup before returning.
	if (!BinBoundarySegments(boundary))
	{
		MemoryCleanup();
		return false;
	}

	if (!GenerateNodes(boundary, outNodes, outBoundaryNodes, rayCastPadding))
	{
		MemoryCleanup();
		return false;
	}

	GenerateElements(outNodes, outRectList);

	MemoryCleanup();

	if (outRectList.size() < 1)
	{
		LogMan::Log("Failed to generate mesh! (Generated zero elements)", LOG_SUCCESS);
		return false;
	}

	LogMan::Log("Succesfully generated mesh!", LOG_SUCCESS);
	return true;
}
