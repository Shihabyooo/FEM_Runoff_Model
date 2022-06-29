#include <string>
#include "DelauneyTriangulation.hpp"


int extVerts[3];//cached external vertices.
int lastID = 0;

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

void GenerateSuperTriangle(std::vector<Vector2D> &nodesList, std::unordered_map<int, Triangle> * outTrianglesList)
{
	LogMan::Log("Generating Super Triangle.");
	auto boundingBox = ComputeBoundingBox(nodesList);

	Vector2D dimensions = boundingBox[1] - boundingBox[0];
	float boundHalfWidth = (dimensions.x / 2.0f) + SUPER_TRIANGLE_PADDING;
	float boundHeight = dimensions.y + 2.0f * SUPER_TRIANGLE_PADDING;

	Vector2D superVert1Pos(	boundingBox[0].x - SUPER_TRIANGLE_PADDING + boundHalfWidth,
							boundingBox[1].y + SUPER_TRIANGLE_PADDING + boundHeight);
	Vector2D superVert2Pos(	boundingBox[0].x - SUPER_TRIANGLE_PADDING - boundHalfWidth,
							boundingBox[0].y - SUPER_TRIANGLE_PADDING);
	Vector2D superVert3Pos(	boundingBox[1].x + SUPER_TRIANGLE_PADDING + boundHalfWidth,
							boundingBox[0].y - SUPER_TRIANGLE_PADDING);

	nodesList.push_back(superVert1Pos);
	nodesList.push_back(superVert2Pos);
	nodesList.push_back(superVert3Pos);

	extVerts[0] = nodesList.size() - 3;
	extVerts[1] = nodesList.size() - 2;
	extVerts[2] = nodesList.size() - 1;

	outTrianglesList->insert({0, Triangle(0, extVerts[0],extVerts[1], extVerts[2], nodesList)});
}

void OptimizeTriangulation(int pivotVertexID, std::vector<Vector2D> const & nodesList, std::vector<Triangle *> & trianglesToOptimize, std::unordered_map<int, Triangle> * trianglesList)
{
	//TODO we don't actually need the shared verts, rather the non shared one of the neighbour tri. Redo related methods to get it directly
	int * sharedVertsIDs = new int[2];

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
					int distantVertID = oldTri->second.GetThirdVertexID(sharedVertsIDs[0], sharedVertsIDs[1]); 

					int newVerts1[3]{ pivotVertexID, distantVertID, testedTri->vertIDs[1] };
					int newVerts2[3]{ pivotVertexID, distantVertID, testedTri->vertIDs[2] };
					Vector2D newNodes1[3]{ nodesList[pivotVertexID], nodesList[distantVertID], nodesList[testedTri->vertIDs[1]]};
					Vector2D newNodes2[3]{ nodesList[pivotVertexID], nodesList[distantVertID], nodesList[testedTri->vertIDs[2]] };

					testedTri->UpdateGeometry(newVerts1, newNodes1);
					oldTri->second.UpdateGeometry(newVerts2, newNodes2);

					//add those new two tris to trianglesToOptimize to test them with distant triangles.
					trianglesToOptimize.push_back(testedTri);
					trianglesToOptimize.push_back(&(oldTri->second));
				}
				break; //stop testing neighbours move on to the next tested triangle.
			}
		}
	}

	delete[] sharedVertsIDs;
}

void DelauneyTriangulation(int vertexID, std::vector<Vector2D> const & nodesList, std::unordered_map<int, Triangle> * trianglesList)
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
	LogMan::Log("Could not fit point " + std::to_string(vertexID) + " in a mesh. Colinear point?", LOG_WARN);
}

std::vector<int> * RemoveExteriorTriangles(std::unordered_map<int, Triangle> * trianglesList) //returns a list of exterior nodes within the mesh
{
	int initialCount = trianglesList->size();
	LogMan::Log("Cleaning exterior triangles.");

	std::vector<int> * outerNodes = new std::vector<int>();
	std::unordered_set<int> tempOuterNodes;

	//loop over triangles, for triangles that test positive for external, add internal nodes (i.e. at mesh edge) to tempOuterNodes, then delete triangle.
	for (auto it = trianglesList->begin(); it != trianglesList->end(); ++it)
	{
		int otherVertices[3];
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
		outerNodes->push_back(*it);

	LogMan::Log("Removed" + std::to_string(initialCount - trianglesList->size()) + " triangles.");

	return outerNodes;
}

bool Triangulate(std::vector<Vector2D> nodesList, std::unordered_map<int, Triangle> * outTrianglesList, std::vector<int> * outBoundaryNodes)
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

	LogMan::Log("Starting triangulation!");
	lastID = 0;

	GenerateSuperTriangle(nodesList, outTrianglesList);

	int nodeCount = nodesList.size() - 3; //actual count, excluding ghost nodes of superTriangle

	for (int i = 0; i < nodeCount; i++)
		DelauneyTriangulation(i, nodesList, outTrianglesList);
	
	outBoundaryNodes = RemoveExteriorTriangles(outTrianglesList);
	
	LogMan::Log("Finished triangulation!", LOG_SUCCESS);
	return true;
}