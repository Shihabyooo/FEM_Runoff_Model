#pragma once
#include "Globals.hpp"

//TODO rewrite this (and delauney triangluation) to have Triangle inherit from Element.
class Triangle
{
public:
	Triangle();
	Triangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, std::vector<Vector2D> const & nodesList);
	Triangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, Vector2D const & node1, Vector2D const & node2, Vector2D const & node3);
	Triangle(size_t _id, size_t const vertexIDs[3], Vector2D const _nodes[3]);
	~Triangle();

	Triangle & operator= (Triangle const & tri2);

	void Subdivide(size_t centroidID, Vector2D centroidPos, size_t baseTriangleID, std::vector<Triangle> & outNewTriangles) const;
	bool ContainsPoint(Vector2D const &point) const;
	
	bool ContainsEdge(size_t vertID1, size_t vertID2) const;
	bool ContainsVertex(size_t vertexID) const; //Contains vertex as a nodal vertix, not inside area. Use ContainsPoint for the latter.

	int GetThirdVertexID(size_t vertID1, size_t vertID2) const; //Warning! Assumes both vertices are indeed nodes of this triangle. Will return first vert doesn't match supplied vertices. Returns -4 if error
	bool IsExternalTriangle(size_t * externalVertices, size_t * outMeshEdgeVerts, int * outMeshEdgeContribCount) const;
	
	bool IsNeighbour(Triangle const &triangle, size_t outsideVertex, size_t * outSharedVertsIDs) const; //sharedVertsIDs = int[2]
	bool IsInsideCircumcircle(Vector2D point) const;
	
	void UpdateGeometry(size_t const vertixIDs[3], Vector2D const _nodes[3]);
	void UpdateGeometry(size_t vertexID1, size_t vertexID2, size_t vertexID3, Vector2D const & node1, Vector2D const & node2, Vector2D const & node3);
	
	double Determinant(Vector2D const &a, Vector2D const &b) const;
	double ComputeArea(Vector2D node1, Vector2D node2, Vector2D node3) const;
	
	Vector2D Centroid() const;

	bool Validate(); //Recomputes area, returns true if area > 0.0, else false.
	void DebugPrintDetails();

	size_t vertIDs[3];
	Vector2D nodes[3]; //TODO It would be much cheaper to store a pointer to (a const)  std::vector<Vector2D> nodes, and reimplement Nodes(int internalVertID) to read the node coord from that vector.
	double area;
	size_t id;
	
	//TODO the value bellow is better stored in a some object in ModelInterface. Leaving it here until I fix the element ID ordering to\
	be sequential starting from 0.
	double elementPrecipitation = 0.0; //in m/hr. The precipitation at last pass.	
};