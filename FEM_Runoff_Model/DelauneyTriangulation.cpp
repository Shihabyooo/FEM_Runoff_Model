#include "DelauneyTriangulation.hpp"

//int extVert1, extVert2, extVert3; //cached external vertices.
int extVerts[3];//cached external vertices.
int lastID = 0;

std::unique_ptr<Vector2[]> ComputeBoundingBox(std::vector<Vector2> &nodesList) //SW then NE corner
{
	auto boundingBox = std::make_unique<Vector2[]>(2);

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

void GenerateSuperTriangle(std::vector<Vector2> &nodesList, std::unordered_map<int, Triangle> * outTrianglesList)
{
	auto boundingBox = ComputeBoundingBox(nodesList);

	Vector2 dimensions = boundingBox[1] - boundingBox[0];
	float boundHalfWidth = (dimensions.x / 2.0f) + SUPER_TRIANGLE_PADDING;
	float boundHeight = dimensions.y + 2.0f * SUPER_TRIANGLE_PADDING;

	Vector2 superVert1Pos(	boundingBox[0].x - SUPER_TRIANGLE_PADDING + boundHalfWidth,
							boundingBox[1].y + SUPER_TRIANGLE_PADDING + boundHeight);
	Vector2 superVert2Pos(	boundingBox[0].x - SUPER_TRIANGLE_PADDING - boundHalfWidth,
							boundingBox[0].y - SUPER_TRIANGLE_PADDING);
	Vector2 superVert3Pos(	boundingBox[1].x + SUPER_TRIANGLE_PADDING + boundHalfWidth,
							boundingBox[0].y - SUPER_TRIANGLE_PADDING);

	nodesList.push_back(superVert1Pos);
	nodesList.push_back(superVert2Pos);
	nodesList.push_back(superVert3Pos);

	extVerts[0] = nodesList.size() - 3;
	extVerts[1] = nodesList.size() - 2;
	extVerts[2] = nodesList.size() - 1;

	outTrianglesList->insert({0, Triangle(0, extVerts[0],extVerts[1], extVerts[2], nodesList)});
	
	std::cout << "\nGenerated supertri\n";
	Print(superVert1Pos);
	Print(superVert2Pos);
	Print(superVert3Pos);
}

void OptimizeTriangulation(int pivotVertexID, std::vector<Vector2> const & nodesList, std::vector<Triangle> & newTriangles, std::unordered_map<int, Triangle> * trianglesList)
{
	int * sharedVertsIDs = NULL;

	for (auto newTri = newTriangles.begin(); newTri != newTriangles.end(); ++newTri)
	{
		for (auto oldTri = trianglesList->begin(); oldTri != trianglesList->end(); ++oldTri)
		{
			if (oldTri->second.IsNeighbour(*newTri, pivotVertexID, &sharedVertsIDs))
			{
				if (oldTri->second.IsInsideCircumcircle(pivotVertexID, nodesList))
				{
					int distantVertID = oldTri->second.GetThirdVertexID(sharedVertsIDs[0], sharedVertsIDs[1]);

					newTri->UpdateGeometry(pivotVertexID, distantVertID, newTri->vertIDs[1], nodesList);
					oldTri->second.UpdateGeometry(pivotVertexID, distantVertID, newTri->vertIDs[2], nodesList);
					
					delete[] sharedVertsIDs;
				}
				break;
			}
		}
	}
}

void DelauneyTriangulation(int vertexID, std::vector<Vector2> const & nodesList, std::unordered_map<int, Triangle> * trianglesList)
{
	//std::cout << "Not implemented yet! " << std::endl; //TODO implement
	std::cout << "------------------------------------------------\n";
	std::cout << "Delauney triangulation for point: " << vertexID << std::endl;

	for (auto it = trianglesList->begin(); it != trianglesList->end(); ++it)
	{
		if (it->second.ContainsPoint(nodesList[vertexID], nodesList))
		{
			std::cout << "Found containing triangle of id: " << it->second.id << std::endl;

			std::vector<Triangle> newTriangles;
			int baseID = lastID++;
			
			//subdivide
			it->second.Subdivide(vertexID, baseID, nodesList, newTriangles);
			
			//remove big triangle
			trianglesList->erase(it->second.id);
			
			//optimize smaller triangles
			OptimizeTriangulation(vertexID, nodesList, newTriangles, trianglesList);

			//add smaller triangles
			for (auto it2 = newTriangles.begin(); it2 != newTriangles.end(); ++it2)
			{
				trianglesList->insert({ baseID, *it2});
				baseID++;
			}
			lastID = baseID;
			return;
		}
	}
}

std::vector<int> * RemoveExteriorTriangles(std::unordered_map<int, Triangle> * trianglesList) //returns a list of exterior nodes within the mesh
{
	std::cout << "Cleaning exterior triangles. Initial count: " << trianglesList->size() << std::endl;
	std::cout << "exterior verts ids: " << extVerts[0] << ", " << extVerts[1] << ", " << extVerts[2] << std::endl;

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

	std::cout << "Count after cleaning: " << trianglesList->size() << std::endl;
	return outerNodes;

}

//void Triangulate(std::vector<Vector2> const &nodesList, std::unordered_map<int, Triangle> * outTrianglesList)
bool Triangulate(std::vector<Vector2> nodesList, std::unordered_map<int, Triangle> * outTrianglesList, std::vector<int> * outBoundaryNodes)
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
		std::cout << "ERROR! Cannot triangulate less than " << MIN_NODES_TO_TRIANGULATE << " nodes." << std::endl;
		return false;
	}

	GenerateSuperTriangle(nodesList, outTrianglesList);

	int nodeCount = nodesList.size() - 3;

	for (int i = 0; i < nodeCount; i++)
	{
		DelauneyTriangulation(i, nodesList, outTrianglesList);

		std::cout << "Triangles so far: " << outTrianglesList->size() << std::endl;
	}

	std::cout << "Finished triangulation!\n\n";
	outBoundaryNodes = RemoveExteriorTriangles(outTrianglesList);
	return true;
}