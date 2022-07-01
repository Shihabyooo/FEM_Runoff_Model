#pragma once
#include <vector>

#include "Globals.hpp"

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

	Vector2D Centroid() const;

	bool Validate(); //Recomputes area, returns true if area > 0.0, else false.
	void DebugPrintDetails();

public:
	int vertIDs[3];
	Vector2D nodes[3];
	double area;
	int id;
	//TODO the values bellow are better stored in a some object in ModelInterface. Leaving them here until I fix the element ID ordering to\
	be sequential starting from 0.
	double elementPrecipitation = 0.0; //in mm. The precipitation at last pass.
	double manningCoef = 0.0; 
	double slopeX = 0.0;
	double slopeY = 0.0;
	double qX = 0.0; //in m3/hr. X-Discharge of last pass;
	double qY = 0.0; //in m3/hr. Y-Discharge of last pass;

private:
	
	double Determinant(Vector2D const &a, Vector2D const &b) const;
	double ComputeArea(Vector2D node1, Vector2D node2, Vector2D node3) const;
};