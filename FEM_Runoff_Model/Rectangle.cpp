#include "Rectangle.hpp"

Rectangle::Rectangle() : Element()
{
	vertCount = 4;
	vertIDs = new size_t[4];
	id = 0;
}

Rectangle::Rectangle(Rectangle const & rect2) : Element(rect2.nodesList)
{
	vertIDs = new size_t[4];
	*this = rect2;
}

Rectangle::Rectangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, size_t vertexID4, std::vector<Vector2D> const * const _nodesList) : Element(_nodesList)
{
	vertCount = 4;
	vertIDs = new size_t[4];
	id = _id;
	size_t vertexIDs[4]{ vertexID1 , vertexID2 , vertexID3 ,vertexID4 };
	UpdateGeometry(vertexIDs);
}

Rectangle::Rectangle(size_t _id, size_t const vertexIDs[4], std::vector<Vector2D> const * const _nodesList) : Element(_nodesList)
{
	vertCount = 4;
	vertIDs = new size_t[4];
	id = _id;
	UpdateGeometry(vertexIDs);
}

Rectangle::~Rectangle()
{
}

Rectangle & Rectangle::operator=(Rectangle const & rect2)
{
	rect = rect2.rect;

	id = rect2.id;
	vertCount = rect2.vertCount; //pointless? vertCount is always 4 for rectangles...
	area = rect2.area;

	if (rect2.vertIDs != NULL)
	{
		for (int i = 0; i < 4; i++)
			vertIDs[i] = rect2.vertIDs[i];
	}
	else if (vertIDs != NULL)
	{
		delete[] vertIDs;
		vertIDs = NULL;
	}
	
	return *this;
}

void Rectangle::UpdateGeometry(size_t const vertixIDs[4]) //TODO This method is ugly. Make it better!
{
	if (vertIDs == NULL)
		vertIDs = new size_t[4];

	vertIDs[0] = vertixIDs[0];
	vertIDs[1] = vertixIDs[1];
	vertIDs[2] = vertixIDs[2];
	vertIDs[3] = vertixIDs[3];

	//divide rect to two tris. Start testing around vert[0], compute three dimensions from vert[0] to other verts.
	//the largest of the three is the hypotenuse (rect's diagonal). The other two are the width and length of rect
	double dim1 = Node(0).DistanceTo(Node(1));
	double dim2 = Node(0).DistanceTo(Node(2));
	double dim3 = Node(0).DistanceTo(Node(3));
	
	//Probably a much better way to do this...
	double height = Min(Min(dim1, dim2), dim3);
	double width = dim1 == height ? (Min(dim2, dim3)) : (dim2 == height? Min (dim1, dim3) : Min(dim1, dim2)); 

	rect = Rect(Vector2D(0, 0), width, height);
	area = rect.Area();

	//create a triangle from vert[0] and the two other verts (other than distantVert), order them in that triangle, then insert distant vert\
	after the second vert.
	int distantVertInternalID = dim1 > dim2 ? (dim1 > dim3 ? 1 : 3) : (dim2 > dim3 ? 2 : 3);
	int otherVert1 = distantVertInternalID == 1 ? 2 : 1;
	int otherVert2 = otherVert1 == 2 ? 3 : (distantVertInternalID == 2 ? 3 : 2);
	
	double subTriArea =  0.5 *	(Node(0).x * (Node(otherVert1).y - Node(otherVert2).y)
		+ Node(otherVert1).x * (Node(otherVert2).y - Node(0).y)
		+ Node(otherVert2).x * (Node(0).y - Node(otherVert1).y));

	//Note: already have vertIDs "backed up" in the recieved VertexIDs, so we can directly assign correct order
	//no change on first vertex
	if (subTriArea < 0.0)
	{
		vertIDs[1] = vertixIDs[otherVert2];
		vertIDs[2] = vertixIDs[distantVertInternalID];
		vertIDs[3] = vertixIDs[otherVert1];
	}
	else
	{
		vertIDs[1] = vertixIDs[otherVert1];
		vertIDs[2] = vertixIDs[distantVertInternalID];
		vertIDs[3] = vertixIDs[otherVert2];
	}
}

double Rectangle::Width()
{
	return rect.Width();
}

double Rectangle::Height()
{
	return rect.Height();
}
