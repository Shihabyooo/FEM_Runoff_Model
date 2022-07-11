#pragma once
#include "Globals.hpp"

//TODO rewrite this (and delauney triangluation) to have Triangle inherit from Element.
class Triangle
{
public:
	Triangle();
	Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2D> const & nodesList);
	Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, Vector2D const & node1, Vector2D const & node2, Vector2D const & node3);
	Triangle(int _id, int const vertexIDs[3], Vector2D const _nodes[3]);
	~Triangle();

	Triangle & operator= (Triangle const & tri2);

	void Subdivide(int centroidID, Vector2D centroidPos, int baseTriangleID, std::vector<Triangle> & outNewTriangles) const;
	bool ContainsPoint(Vector2D const &point) const;
	
	bool ContainsEdge(int vertID1, int vertID2) const;
	bool ContainsVertex(int vertexID) const; //Contains vertex as a nodal vertix, not inside area. Use ContainsPoint for the latter.

	int GetThirdVertexID(int vertID1, int vertID2) const; //Warning! Assumes both vertices are indeed nodes of this triangle. Will return first vert doesn't match supplied vertices. Returns -4 if error
	bool IsExternalTriangle(int * externalVertices, int * outMeshEdgeVerts, int * outMeshEdgeContribCount) const;
	
	bool IsNeighbour(Triangle const &triangle, int outsideVertex, int * outSharedVertsIDs) const; //sharedVertsIDs = int[2]
	bool IsInsideCircumcircle(Vector2D point) const;
	
	void UpdateGeometry(int const vertixIDs[3], Vector2D const _nodes[3]);
	void UpdateGeometry(int vertexID1, int vertexID2, int vertexID3, Vector2D const & node1, Vector2D const & node2, Vector2D const & node3);
	
	double Determinant(Vector2D const &a, Vector2D const &b) const;
	double ComputeArea(Vector2D node1, Vector2D node2, Vector2D node3) const;
	
	Vector2D Centroid() const;

	bool Validate(); //Recomputes area, returns true if area > 0.0, else false.
	void DebugPrintDetails();

	int vertIDs[3];
	Vector2D nodes[3]; //TODO It would be much cheaper to store a pointer to (a const)  std::vector<Vector2D> nodes, and reimplement Nodes(int internalVertID) to read the node coord from that vector.
	double area;
	int id;
	
	//TODO the value bellow is better stored in a some object in ModelInterface. Leaving it here until I fix the element ID ordering to\
	be sequential starting from 0.
	double elementPrecipitation = 0.0; //in m/hr. The precipitation at last pass.	
};