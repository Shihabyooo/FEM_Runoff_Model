#pragma once
#include <vector>
#include <MatricesPP.hpp>

//#include <iostream>

#include "Globals.hpp"

class Triangle
{
public:
	Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2> const & nodesList);
	~Triangle();

	void Subdivide(int centroidID, int baseTriangleID, std::vector<Vector2> const & nodesList, std::vector<Triangle> & outNewTriangles);
	bool ContainsPoint(Vector2 const &point, std::vector<Vector2> const & nodesList);
	
	bool ContainsEdge(int vertID1, int vertID2) const;
	bool ContainsVertex(int vertexID) const;

	int GetThirdVertexID(int vertID1, int vertID2) const; //Warning! Assumes both vertices are indeed nodes of this triangle. Will return first vert doesn't match supplied vertices. Returns -4 if error
	bool IsExternalTriangle(int * externalVertices, int * outMeshEdgeVerts, int * outMeshEdgeContribCount);
	
	bool IsNeighbour(Triangle const &triangle, int outsideVertex, int * outSharedVertsIDs); //sharedVertsIDs = int[2]
	bool IsInsideCircumcircle(int vertexID, std::vector<Vector2> const & nodesList);
	
	void UpdateGeometry(int const * vertixIDs, std::vector<Vector2> const & _nodesList);
	void UpdateGeometry(int const & vertexID1, int const & vertexID2, int const & vertexID3, std::vector<Vector2> const & _nodesList);

	Vector2 Node(int internalVertexID, std::vector<Vector2> const & nodesList) const;

	void DebugPrintDetails();

public:
	int vertIDs[3];
	double area;
	int id;

private:
	

	double Determinant(Vector2 const &a, Vector2 const &b);
	double ComputeArea(Vector2 node1, Vector2 node2, Vector2 node3);

//private:
//	static std::vector<Vector2> * nodesList ;
};