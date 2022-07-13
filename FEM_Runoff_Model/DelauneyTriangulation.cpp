#include <string>
#include <chrono>
#include "DelauneyTriangulation.hpp"

//int extVerts[3];//cached external vertices.
size_t extVerts[4];//cached external vertices.
int lastID = 0;
std::vector<int> failedPoints; //to hold points that failed to triangulate on first pass.

std::unique_ptr<Vector2D[]> ComputeBoundingBox(std::vector<Vector2D> &nodesList) //SW then NE corner
{	//TODO replace this with a global func (since it's also computed in ModelInterface and required in several places).
	auto boundingBox = std::make_unique<Vector2D[]>(2);

	boundingBox[0] = nodesList[0];
	boundingBox[1] = nodesList[0];

	for (auto it = nodesList.begin(); it < nodesList.end(); it++)
	{
		boundingBox[0].x = Min(boundingBox[0].x, it->x);
		boundingBox[0].y = Min(boundingBox[0].y, it->y);

		boundingBox[1].x = Max(boundingBox[1].x, it->x);
		boundingBox[1].y = Max(boundingBox[1].y, it->y);
	}

	return boundingBox;
}

void GenerateSuperTriangle(std::vector<Vector2D> &nodesList, double padding, std::unordered_map<size_t, Triangle> * outTrianglesList, Triangle * outSuperTriangles)
{
	LogMan::Log("Generating Super Triangle.");
	auto boundingBox = ComputeBoundingBox(nodesList);

	/*Vector2D dimensions = boundingBox[1] - boundingBox[0];
	double padding = SUPER_TRIANGLE_PADDING * Max(dimensions.x, dimensions.y);*/
	/*float boundHalfWidth = (dimensions.x / 2.0f) + SUPER_TRIANGLE_PADDING;
	float boundHeight = dimensions.y + 2.0f * SUPER_TRIANGLE_PADDING;*/

	Vector2D superVert1Pos(	boundingBox[1].x + padding,
							boundingBox[0].y - padding);
	Vector2D superVert2Pos(	boundingBox[0].x - padding,
							boundingBox[1].y + padding);
	Vector2D superVert3Pos(	boundingBox[0].x - padding,
							boundingBox[0].y - padding);
	Vector2D superVert4Pos(	boundingBox[1].x + padding,
							boundingBox[1].y + padding);

	nodesList.push_back(superVert1Pos);
	nodesList.push_back(superVert2Pos);
	nodesList.push_back(superVert3Pos);
	nodesList.push_back(superVert4Pos);

	extVerts[0] = nodesList.size() - 4;
	extVerts[1] = nodesList.size() - 3;
	extVerts[2] = nodesList.size() - 2;
	extVerts[3] = nodesList.size() - 1;

	outTrianglesList->insert({-1, Triangle(-1, extVerts[0],extVerts[1], extVerts[2], nodesList)});
	outTrianglesList->insert({-2, Triangle(-2, extVerts[0],extVerts[1], extVerts[3], nodesList)});

	outSuperTriangles[0] = Triangle(-1, extVerts[0], extVerts[1], extVerts[2], nodesList);
	outSuperTriangles[1] = Triangle(-2, extVerts[0], extVerts[1], extVerts[3], nodesList);
}

void OptimizeTriangulation(int pivotVertexID, std::vector<Vector2D> const & nodesList, std::vector<Triangle *> & trianglesToOptimize, std::unordered_map<size_t, Triangle> * trianglesList)
{
	//TODO we don't actually need the shared verts, rather the non shared one of the neighbour tri. Redo related methods to get it directly
	size_t * sharedVertsIDs = new size_t[2];

	while (trianglesToOptimize.size() > 1)
	{
		Triangle * testedTri = trianglesToOptimize.back();
		trianglesToOptimize.pop_back();

		//search for neighbour
		for (auto oldTri = trianglesList->begin(); oldTri != trianglesList->end(); ++oldTri)
		{
			if (testedTri->id != oldTri->second.id && testedTri->IsNeighbour(oldTri->second, pivotVertexID, sharedVertsIDs)) //this is a neighbour
			{
				if (oldTri->second.IsInsideCircumcircle(nodesList[pivotVertexID])) //This is a viable optimization candidate.
				{
					size_t distantVertID = oldTri->second.GetThirdVertexID(sharedVertsIDs[0], sharedVertsIDs[1]);

					size_t newVerts1[3]{ pivotVertexID, distantVertID, sharedVertsIDs[0] };
					size_t newVerts2[3]{ pivotVertexID, distantVertID, sharedVertsIDs[1] };

					Vector2D newNodes1[3]{ nodesList[pivotVertexID], nodesList[distantVertID], nodesList[sharedVertsIDs[0]] };
					Vector2D newNodes2[3]{ nodesList[pivotVertexID], nodesList[distantVertID], nodesList[sharedVertsIDs[1]] };

					/*std::cout << "tested tris:\n";
					testedTri->DebugPrintDetails();
					oldTri->second.DebugPrintDetails();
					std::cout << "pivot vert: " << pivotVertexID << std::endl;
					std::cout << "distant vert: " << distantVertID << std::endl;
					std::cout << "edgeverts: " << sharedVertsIDs[0] << ", " << sharedVertsIDs[1] << std::endl;
					std::cout << "new Verts1: " << newVerts1[0] << ", " << newVerts1[1] << ", " << newVerts1[2] << std::endl;
					std::cout << "new Verts2: " << newVerts2[0] << ", " << newVerts2[1] << ", " << newVerts2[2] << std::endl;*/

					testedTri->UpdateGeometry(newVerts1, newNodes1);
					oldTri->second.UpdateGeometry(newVerts2, newNodes2);

					//add those new two tris to trianglesToOptimize to test them with distant triangles.
					trianglesToOptimize.push_back(testedTri);
					trianglesToOptimize.push_back(&(oldTri->second));

					////test
					//if (!testedTri->Validate())
					//{
					//	std::cout << "!!!!!!!!!!!!!!!!!!Error at:";
					//	testedTri->DebugPrintDetails();
					//}
					//if (!oldTri->second.Validate())
					//{
					//	std::cout << "!!!!!!!!!!!!!!!!!!Error at:";
					//	oldTri->second.DebugPrintDetails();
					//}
					////endtest
				}
				break; //stop testing neighbours move on to the next tested triangle.
			}
		}
	}

	delete[] sharedVertsIDs;
}

void DelauneyTriangulation(int vertexID, std::vector<Vector2D> const & nodesList, std::unordered_map<size_t, Triangle> * trianglesList, bool firstPass = true)
{
	for (auto it = trianglesList->begin(); it != trianglesList->end(); ++it)
	{
		if (it->second.ContainsPoint(nodesList[vertexID]))
		{
			int baseID = lastID++;
			std::vector<Triangle> newTriangles;
			
			//subdivide
			it->second.Subdivide(vertexID, nodesList[vertexID], baseID, newTriangles);

			//remove big triangle
			trianglesList->erase(it->second.id);
			
			//optimize smaller triangles
			//create a list of pointers to triangels that include current vertexID. For now it's exactly the contents
			//of newTriangles, but more will be added inside OptimizeTriangulation(), so we keep our newTriangles list separate.
			std::vector<Triangle *> trisToOptimize;
			for (auto it2 = newTriangles.begin(); it2 != newTriangles.end(); ++it2)
				trisToOptimize.push_back(&(*it2));

			OptimizeTriangulation(vertexID, nodesList, trisToOptimize, trianglesList);

			//add smaller triangles to map
			for (auto it2 = newTriangles.begin(); it2 != newTriangles.end(); ++it2)
			{
				if (!trianglesList->insert({ it2->id, *it2 }).second)
					LogMan::Log("Failed to insert triangle with ID: " + it2->id, LOG_ERROR);
			}
			
			lastID = baseID + 3;
			return;
		}
	}
	
	if (firstPass)
		failedPoints.push_back(vertexID);
	else
		LogMan::Log("Could not fit point " + std::to_string(vertexID) + " in the mesh.", LOG_WARN); 
}

void OptimizeAll(std::vector<Vector2D> const & nodesList, std::unordered_map<size_t, Triangle> * triangles)
{
	LogMan::Log("Global Optimization run");

	for (int i = 0; i < nodesList.size() - 4; i++)
	{
		std::vector<Triangle *> toOptimize;
		for (auto it = triangles->begin(); it != triangles->end(); ++it)
		{
			if (it->second.ContainsVertex(i))
				toOptimize.push_back(&(it->second));
		}
		OptimizeTriangulation(i, nodesList, toOptimize, triangles);
	}
}

void RemoveExteriorTriangles(std::unordered_map<size_t, Triangle> * trianglesList, std::vector<size_t> * outBoundaryNodes) //returns a list of exterior nodes within the mesh
{
	int initialCount = trianglesList->size();
	LogMan::Log("Cleaning exterior triangles.");

	//std::vector<int> * outerNodes = new std::vector<int>();
	std::unordered_set<size_t> tempOuterNodes;

	//loop over triangles, for triangles that test positive for external, add internal nodes (i.e. at mesh edge) to tempOuterNodes, then delete triangle.
	for (auto it = trianglesList->begin(); it != trianglesList->end(); ++it)
	{
		size_t otherVertices[3];
		int meshEdgeVertContrib;
		if (it->second.IsExternalTriangle(extVerts, otherVertices, &meshEdgeVertContrib))
		{
			for (int i = 0; i < meshEdgeVertContrib; i++)
				tempOuterNodes.insert(otherVertices[i]);
			
			trianglesList->erase(it--);
		}
	}

	//copy tempOuterNodes to outerNodes
	for (auto it = tempOuterNodes.begin(); it != tempOuterNodes.end(); ++it)
		outBoundaryNodes->push_back(*it);

	LogMan::Log("Removed" + std::to_string(initialCount - trianglesList->size()) + " triangles.");
}

void Cleanup(std::vector<Vector2D> const & originalNodesList, std::unordered_map<size_t, Triangle> * trianglesList)
{
	//Two things to do: Remove invalid tris, and have ID sequentially from 0.

	//Because of the carried out above, we may have created triangles that would, when referenced to the original node, have\
	an area of zero (because, e.g. all three points are colinear.

	//Note: while the resulting IDs would be sequential, they would not necessarily be ordered in a way that makes sense.

	LogMan::Log("Internal cleanup pass");
	size_t counter = 0;
	std::unordered_map<size_t, Triangle> cleanedTriList;

	for (auto it = trianglesList->begin(); it != trianglesList->end(); ++it)
	{
		size_t const * verts = it->second.vertIDs;
		Vector2D originalNodes[3]{ originalNodesList[verts[0]], originalNodesList[verts[1]], originalNodesList[verts[2]] };
		it->second.UpdateGeometry(verts, originalNodes);
		
		if (!it->second.Validate()) //remove
		{
			LogMan::Log("Removed invalid triangle: " + std::to_string(it->second.id));
		}
		else //update id and copy to cleaned list
		{
			it->second.id = counter; //Adjust in place before copying to cleanedTriList. Ok since old data will be deleted anyway.
			cleanedTriList.insert({ counter, it->second });
			counter++;
		}
	}
	
	*trianglesList = std::move(cleanedTriList);
}

bool Triangulate(std::vector<Vector2D> const & nodesList, double superTrianglePadding, std::unordered_map<size_t, Triangle> * outTrianglesList, std::vector<size_t> * outBoundaryNodes, Triangle * outSuperTriangles)
{
	//check if nodesList has enough elements, else return false.
	//compute bounding box
	//Generate super triangle
	//insert verts of super triangle at the end of the nodesList (so order of original verts won't change)
	//loop over each vert in nodesList (except last three)
		//Call DelauneyTriangulation() for each node.
	//Call RemoveExteriorTriangles()
	//Compute edge nodes
	//Set and return computation results

	if (nodesList.size() < MIN_NODES_TO_TRIANGULATE)
	{
		LogMan::Log(("Cannot triangulate less than than: " + std::to_string(MIN_NODES_TO_TRIANGULATE) + " nodes"), LOG_ERROR);
		return false;
	}

	std::vector<Vector2D> tempNodes = nodesList; //make a copy of the nodes that we can modify

	LogMan::Log("Starting triangulation!");
	lastID = 0;
	failedPoints.clear();
	GenerateSuperTriangle(tempNodes, superTrianglePadding, outTrianglesList, outSuperTriangles);

	int nodeCount = tempNodes.size() - 4; //actual count, excluding ghost nodes of superTriangle

	//first pass
	LogMan::Log("First pass");
	for (int i = 0; i < nodeCount; i++)
	{
		//LogMan::Log("At node: " + std::to_string(i));
		DelauneyTriangulation(i, tempNodes, outTrianglesList);
	}

	//second pass	
	LogMan::Log("Second pass");
	//loop over failed points
	srand(time(0));
	for (auto it = failedPoints.begin(); it != failedPoints.end(); ++it)
	{
		//add some jitter to ths point
		tempNodes[*it].x += tempNodes[*it].x * static_cast<double>(rand() % 10) * 0.00001 * pow(-1.0, rand() % 2);
		srand(tempNodes[*it].y);
		tempNodes[*it].y += tempNodes[*it].y * static_cast<double>(rand() % 10) * 0.00001 * pow(-1.0, rand() % 2);
		//triangulate
		DelauneyTriangulation(*it, tempNodes, outTrianglesList, false);
	}
	
	OptimizeAll(tempNodes, outTrianglesList);
	RemoveExteriorTriangles(outTrianglesList, outBoundaryNodes);
	Cleanup(nodesList, outTrianglesList);

	LogMan::Log("Finished triangulation!", LOG_SUCCESS);
	return true;
}