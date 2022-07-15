#pragma once
#include "FEM_Element.hpp"

//TODO rewrite this (and delauney triangluation) to have Triangle inherit from Element.
class Triangle : public Element
{
public:
	Triangle();
	Triangle(Triangle const & tri2);
	Triangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, std::vector<Vector2D> const * const _nodesList);
	Triangle(size_t _id, size_t const vertexIDs[3], std::vector<Vector2D> const * const _nodesList);
	~Triangle();

	Triangle & operator= (Triangle const & tri2);

	void Subdivide(size_t centroidID, size_t baseTriangleID, std::vector<Triangle> & outNewTriangles) const;
	bool ContainsPoint(Vector2D const &point) const override;

	int GetThirdVertexID(size_t vertID1, size_t vertID2) const; //Warning! Assumes both vertices are indeed nodes of this triangle. Will return first vert doesn't match supplied vertices. Returns -4 if error
	bool IsExternalTriangle(size_t * externalVertices, size_t * outMeshEdgeVerts, int * outMeshEdgeContribCount) const;
	
	bool IsNeighbour(Triangle const &triangle, size_t outsideVertex, size_t * outSharedVertsIDs) const; //sharedVertsIDs = int[2]
	bool IsInsideCircumcircle(size_t vertID) const;
	
	void UpdateGeometry(size_t const vertexIDs[3]);
	
	double Determinant(Vector2D const &a, Vector2D const &b) const;
	double ComputeArea() const override;
};