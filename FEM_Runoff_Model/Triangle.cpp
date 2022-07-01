#include "Triangle.hpp"
#include <MatricesPP.hpp>

Triangle::Triangle()
{
}

Triangle::Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2D> const & nodesList)
{
	id = _id;
	UpdateGeometry(vertexID1, vertexID2, vertexID3, nodesList[vertexID1], nodesList[vertexID2], nodesList[vertexID3]);
}

Triangle::Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, Vector2D const & node1, Vector2D const & node2, Vector2D const & node3)
{
	id = _id;
	UpdateGeometry(vertexID1, vertexID2, vertexID3, node1, node2, node3);
}

Triangle::Triangle(int _id, int const vertexIDs[3], Vector2D const _nodes[3])
{
	id = _id;
	UpdateGeometry(vertexIDs, _nodes);
}

Triangle::~Triangle()
{
}

Triangle & Triangle::operator=(Triangle const & tri2)
{
	id = tri2.id;
	area = tri2.area;
	vertIDs[0] = tri2.vertIDs[0];
	vertIDs[1] = tri2.vertIDs[1];
	vertIDs[2] = tri2.vertIDs[2];

	nodes[0] = tri2.nodes[0];
	nodes[1] = tri2.nodes[1];
	nodes[2] = tri2.nodes[2];

	return *this;
}

void Triangle::Subdivide(int centroidID, Vector2D centroidPos, int baseTriangleID, std::vector<Triangle> & outNewTriangles) const
{
	outNewTriangles.push_back(Triangle(baseTriangleID	 , centroidID, vertIDs[0], vertIDs[1], centroidPos, nodes[0], nodes[1]));
	outNewTriangles.push_back(Triangle(baseTriangleID + 1, centroidID, vertIDs[0], vertIDs[2], centroidPos, nodes[0], nodes[2]));
	outNewTriangles.push_back(Triangle(baseTriangleID + 2, centroidID, vertIDs[1], vertIDs[2], centroidPos, nodes[1], nodes[2]));
}

bool Triangle::ContainsPoint(Vector2D const & point) const
{
	//https://mathworld.wolfram.com/TriangleInterior.html

	Vector2D v1 = nodes[1] - nodes[0];
	Vector2D v2 = nodes[2] - nodes[0];

	double a = (Determinant(point, v2) - Determinant(nodes[0], v2)) / Determinant(v1, v2);
	double b = -1.0 * (Determinant(point, v1) - Determinant(nodes[0], v1)) / Determinant(v1, v2);

	if (a > 0.0 && b > 0.0 && (a + b) < 1.0)
		return true;

	return false;
}

bool Triangle::ContainsEdge(int vertID1, int vertID2) const
{
	bool check1 = false, check2 = false;

	for (int i = 0; i < 3; i++)
	{
		check1 = check1 || (vertIDs[i] == vertID1);
		check2 = check2 || (vertIDs[i] == vertID2);
	}

	return (check1 && check2);
}

bool Triangle::ContainsVertex(int vertexID) const
{
	return ((vertIDs[0] == vertexID) || (vertIDs[1] == vertexID) || (vertIDs[2] == vertexID));
}

int Triangle::GetThirdVertexID(int vertID1, int vertID2) const
{
	for (int i = 0; i < 3; i++)
	{
		if (vertIDs[i] != vertID1 && vertIDs[i] != vertID2)
			return vertIDs[i];
	}

	return -4;
}

bool Triangle::IsExternalTriangle(int * externalVertices, int * outMeshEdgeVerts, int * outMeshEdgeContribCount) const
{
	*outMeshEdgeContribCount = 0;
	
	bool isInternal = true;

	for (int i = 0; i < 4; i++)
		isInternal = isInternal && !ContainsVertex(externalVertices[i]);
	

	if (!isInternal)
	{
		for (int i = 0; i < 3; i++)
		{
			if (vertIDs[i] != externalVertices[0] && vertIDs[i] != externalVertices[1] && vertIDs[i] != externalVertices[2])
			{
				outMeshEdgeVerts[*outMeshEdgeContribCount] = vertIDs[i];
				(*outMeshEdgeContribCount)++;
			}
		}
	}
	
	return !isInternal;
}

bool Triangle::IsNeighbour(Triangle const & triangle, int outsideVertex, int * outSharedVertsIDs) const
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

bool Triangle::IsInsideCircumcircle(Vector2D pos) const
{
	//https://stackoverflow.com/questions/39984709/how-can-i-check-wether-a-point-is-inside-the-circumcircle-of-3-points

	Matrix_f32 matrix(3, 3);
	//TODO a, b, and c are remnants from a previou implementationi. Remove them and directly access nodes[] in the lines bellow.
	Vector2D a = nodes[0];
	Vector2D b = nodes[1];
	Vector2D c = nodes[2];

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

void Triangle::UpdateGeometry(int const vertixIDs[3], Vector2D const _nodes[3])
{
	UpdateGeometry(vertixIDs[0], vertixIDs[1], vertixIDs[2], _nodes[0], _nodes[1], _nodes[2]);
}

void Triangle::UpdateGeometry(int vertexID1, int vertexID2, int vertexID3,
								Vector2D const & node1, Vector2D const & node2, Vector2D const & node3)
{
	vertIDs[0] = vertexID1;
	vertIDs[1] = vertexID2;
	vertIDs[2] = vertexID3;

	nodes[0] = node1;
	nodes[1] = node2;
	nodes[2] = node3;

	area = ComputeArea(nodes[0], nodes[1], nodes[2]);
	if (area < 0.0f) //points were passed clockwise, flip than and absolute the area.
	{
		int tempVertID = vertIDs[1];		
		vertIDs[1] = vertIDs[2];
		vertIDs[2] = tempVertID;
		
		Vector2D tempNode = nodes[1];
		nodes[1] = nodes[2];
		nodes[2] = tempNode;

		area = abs(area); //no need to compute area again, vertex reording only changes sign.
	}
}

Vector2D Triangle::Centroid() const
{
	Vector2D centroid;

	centroid = nodes[0] + nodes[1] + nodes[2];
	centroid = Vector2D(centroid.x / 3.0, centroid.y / 3.0);

	return centroid;
}

bool Triangle::Validate()
{
	area = ComputeArea(nodes[0], nodes[1], nodes[2]);
	return area > 0.0;
}

void Triangle::DebugPrintDetails()
{
	std::cout << "-ID: " << id << ".  \tVerts: " << vertIDs[0] << ", " << vertIDs[1] << ", " << vertIDs[2] << " - nodes: ";
	Print(nodes[0], true);
	std::cout << " | ";
	Print(nodes[1], true);
	std::cout << " | ";
	Print(nodes[1], true);
	std::cout << std::endl;
}

double Triangle::Determinant(Vector2D const & a, Vector2D const & b) const
{
	return ((a.x * b.y) - (a.y * b.x));
}

double Triangle::ComputeArea(Vector2D node1, Vector2D node2, Vector2D node3) const
{
	return	0.5 *	( node1.x * (node2.y - node3.y)
					+ node2.x * (node3.y - node1.y)
					+ node3.x * (node1.y - node2.y));
}