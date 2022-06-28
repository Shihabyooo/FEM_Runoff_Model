#pragma once
#include <vector>

#include "Globals.hpp"

class Triangle
{
public:
	Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2> const & nodesList);
	~Triangle();

	void Subdivide(int centroidID, int baseTriangleID, std::vector<Vector2> const & nodesList, std::vector<Triangle> & outNewTriangles) const;
	bool ContainsPoint(Vector2 const &point, std::vector<Vector2> const & nodesList) const;
	
	bool ContainsEdge(int vertID1, int vertID2) const;
	bool ContainsVertex(int vertexID) const;

	int GetThirdVertexID(int vertID1, int vertID2) const; //Warning! Assumes both vertices are indeed nodes of this triangle. Will return first vert doesn't match supplied vertices. Returns -4 if error
	bool IsExternalTriangle(int * externalVertices, int * outMeshEdgeVerts, int * outMeshEdgeContribCount) const;
	
	bool IsNeighbour(Triangle const &triangle, int outsideVertex, int * outSharedVertsIDs) const; //sharedVertsIDs = int[2]
	bool IsInsideCircumcircle(int vertexID, std::vector<Vector2> const & nodesList) const;
	
	void UpdateGeometry(int const * vertixIDs, std::vector<Vector2> const & _nodesList);
	void UpdateGeometry(int const & vertexID1, int const & vertexID2, int const & vertexID3, std::vector<Vector2> const & _nodesList);

	Vector2 Node(int internalVertexID, std::vector<Vector2> const & nodesList) const;

	void DebugPrintDetails();

public:
	int vertIDs[3];
	double area;
	int id;
	//TODO the values bellow are better stored in a some object in ModelInterface. Leaving them here until I fix the element ID ordering to\
	be sequential starting from 0.
	double elementPrecipitation = 0.0; //in mm. The precipitation at last pass.
	double manningCoef = 0.0; 

private:
	
	double Determinant(Vector2 const &a, Vector2 const &b) const;
	double ComputeArea(Vector2 node1, Vector2 node2, Vector2 node3) const;
};