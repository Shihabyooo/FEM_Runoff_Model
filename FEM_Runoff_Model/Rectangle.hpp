#pragma once
#include "FEM_Element.hpp"

class Rectangle : public Element
{
public:
	Rectangle();
	Rectangle(Rectangle const & rect2);
	Rectangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, size_t vertexID4, std::vector<Vector2D> const * const _nodesList);
	Rectangle(size_t _id, size_t const vertexIDs[4], std::vector<Vector2D> const * const _nodesList);
	~Rectangle();

	Rectangle & operator= (Rectangle const & rect2);

	void UpdateGeometry(size_t const vertixIDs[4]);
	
	double ComputeArea() const override;

	double Width();
	double Height();

private:
	Rect rect; //rect is not exactly of the position of this element, because Rect assumes orientation matching x-y axes, while element may be rotated.
};