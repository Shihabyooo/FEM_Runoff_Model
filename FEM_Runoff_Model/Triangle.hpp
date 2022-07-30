#pragma once
#include "FEM_Element.hpp"

class Triangle : public Element
{
public:
	Triangle();
	Triangle(Triangle const & tri2);
	Triangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, std::vector<Vector2D> const * const _nodesList);
	Triangle(size_t _id, size_t const vertexIDs[3], std::vector<Vector2D> const * const _nodesList);
	~Triangle();

	Triangle & operator= (Triangle const & tri2);

	bool ContainsPoint(Vector2D const &point) const override;
	void UpdateGeometry(size_t const vertexIDs[3]);
	double Determinant(Vector2D const &a, Vector2D const &b) const;
	double ComputeArea() const override;
};