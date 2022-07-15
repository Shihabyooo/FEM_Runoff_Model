#include "Triangle.hpp"
//#include <MatricesPP.hpp>

Triangle::Triangle()
{
}

Triangle::Triangle(Triangle const & tri2) : Element(tri2.nodesList)
{
	vertIDs = new size_t[3];
	*this = tri2;
}

Triangle::Triangle(size_t _id, size_t vertexID1, size_t vertexID2, size_t vertexID3, std::vector<Vector2D> const * const _nodesList) : Element(_nodesList)
{
	vertCount = 3;
	vertIDs = new size_t[3];
	id = _id;
	size_t vertexIDs[3]{ vertexID1 , vertexID2 , vertexID3 };
	UpdateGeometry(vertexIDs);
}

Triangle::Triangle(size_t _id, size_t const vertexIDs[3], std::vector<Vector2D> const * const _nodesList) : Element(_nodesList)
{
	vertCount = 3;
	vertIDs = new size_t[3];
	id = _id;
	UpdateGeometry(vertexIDs);
}


Triangle::~Triangle()
{
}

Triangle & Triangle::operator=(Triangle const & tri2)
{
	id = tri2.id;
	vertCount = tri2.vertCount; //pointless? 
	area = tri2.area;

	if (tri2.vertIDs != NULL)
	{
		if (vertIDs == NULL)
			vertIDs = new size_t[3];

		for (int i = 0; i < 3; i++)
			vertIDs[i] = tri2.vertIDs[i];
	}
	else if (vertIDs != NULL) //If the assigned vertIDs is null and this one isn't, delete this one so it can also be null.
	{
		delete[] vertIDs;
		vertIDs = NULL;
	}

	return *this;
}

void Triangle::Subdivide(size_t centroidID, size_t baseTriangleID, std::vector<Triangle> & outNewTriangles) const
{
	outNewTriangles.push_back(Triangle(baseTriangleID	 , centroidID, vertIDs[0], vertIDs[1], nodesList));
	outNewTriangles.push_back(Triangle(baseTriangleID + 1, centroidID, vertIDs[0], vertIDs[2], nodesList));
	outNewTriangles.push_back(Triangle(baseTriangleID + 2, centroidID, vertIDs[1], vertIDs[2], nodesList));
}

bool Triangle::ContainsPoint(Vector2D const & point) const
{
	//https://mathworld.wolfram.com/TriangleInterior.html

	Vector2D v1 = Node(1) - Node(0);
	Vector2D v2 = Node(2) - Node(0);

	double a = (Determinant(point, v2) - Determinant(Node(0), v2)) / Determinant(v1, v2);
	double b = -1.0 * (Determinant(point, v1) - Determinant(Node(0), v1)) / Determinant(v1, v2);

	if (a > 0.0 && b > 0.0 && (a + b) < 1.0)
		return true;

	return false;
}

int Triangle::GetThirdVertexID(size_t vertID1, size_t vertID2) const
{
	for (int i = 0; i < 3; i++)
	{
		if (vertIDs[i] != vertID1 && vertIDs[i] != vertID2)
			return vertIDs[i];
	}

	return -4;
}

bool Triangle::IsExternalTriangle(size_t * externalVertices, size_t * outMeshEdgeVerts, int * outMeshEdgeContribCount) const
{
	*outMeshEdgeContribCount = 0;
	
	bool isInternal = true;

	for (int i = 0; i < 4; i++)
		isInternal = isInternal && !ContainsVertex(externalVertices[i]);
	

	if (!isInternal)
	{
		for (int i = 0; i < 3; i++)
		{
			if (vertIDs[i] != externalVertices[0] && vertIDs[i] != externalVertices[1] && vertIDs[i] != externalVertices[2] && vertIDs[i] != externalVertices[3])
			{
				outMeshEdgeVerts[*outMeshEdgeContribCount] = vertIDs[i];
				(*outMeshEdgeContribCount)++;
			}
		}
	}
	
	return !isInternal;
}

bool Triangle::IsNeighbour(Triangle const & triangle, size_t outsideVertex, size_t * outSharedVertsIDs) const
{
	std::vector<int> edgeVerts;

	for (int i = 0; i < 3; i++)
	{
		if (vertIDs[i] != outsideVertex)
			edgeVerts.push_back(vertIDs[i]);
	}

	if (triangle.ContainsEdge(edgeVerts[0], edgeVerts[1]))
	{
		outSharedVertsIDs[0] = edgeVerts[0];
		outSharedVertsIDs[1] = edgeVerts[1];
		return true;
	}
	else
	{
		return false;
	}
}

bool Triangle::IsInsideCircumcircle(size_t vertID) const
{
	//https://stackoverflow.com/questions/39984709/how-can-i-check-wether-a-point-is-inside-the-circumcircle-of-3-points

	Matrix_f32 matrix(3, 3);

	Vector2D const & a = Node(0);
	Vector2D const & b = Node(1);
	Vector2D const & c = Node(2);
	Vector2D const & pos = (*nodesList)[vertID];

	matrix[0][0] = a.x - pos.x;
	matrix[0][1] = a.y - pos.y;
	matrix[0][2] = pow(a.x - pos.x, 2.0f) + pow(a.y - pos.y, 2.0f);

	matrix[1][0] = b.x - pos.x;
	matrix[1][1] = b.y - pos.y;
	matrix[1][2] = pow(b.x - pos.x, 2.0f) + pow(b.y - pos.y, 2.0f);

	matrix[2][0] = c.x - pos.x;
	matrix[2][1] = c.y - pos.y;
	matrix[2][2] = pow(c.x - pos.x, 2.0f) + pow(c.y - pos.y, 2.0f);

	return (matrix.Determinant() > 0.0f);
}

void Triangle::UpdateGeometry(size_t const vertexIDs[3])
{
	vertIDs[0] = vertexIDs[0];
	vertIDs[1] = vertexIDs[1];
	vertIDs[2] = vertexIDs[2];

	area = ComputeArea();
	if (area < 0.0f) //points were passed clockwise, flip than and absolute the area.
	{
		int tempVertID = vertIDs[1];		
		vertIDs[1] = vertIDs[2];
		vertIDs[2] = tempVertID;
		
		area = abs(area); //no need to compute area again, vertex reording only changes sign.
	}
}

double Triangle::Determinant(Vector2D const & a, Vector2D const & b) const
{
	return ((a.x * b.y) - (a.y * b.x));
}

double Triangle::ComputeArea() const
{
	return	0.5 *	( Node(0).x * (Node(1).y - Node(2).y)
					+ Node(1).x * (Node(2).y - Node(0).y)
					+ Node(2).x * (Node(0).y - Node(1).y));
}