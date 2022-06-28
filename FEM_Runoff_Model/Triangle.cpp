#include "Triangle.hpp"
#include <MatricesPP.hpp>

Triangle::Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2> const & nodesList)
{
	id = _id;
	UpdateGeometry(vertexID1, vertexID2, vertexID3, nodesList);
}

Triangle::~Triangle()
{
}

void Triangle::Subdivide(int centroidID, int baseTriangleID, std::vector<Vector2> const & nodesList, std::vector<Triangle>& outNewTriangles) const
{
	outNewTriangles.push_back(Triangle(baseTriangleID	 , centroidID, vertIDs[0], vertIDs[1], nodesList));
	outNewTriangles.push_back(Triangle(baseTriangleID + 1, centroidID, vertIDs[0], vertIDs[2], nodesList));
	outNewTriangles.push_back(Triangle(baseTriangleID + 2, centroidID, vertIDs[1], vertIDs[2], nodesList));
}

bool Triangle::ContainsPoint(Vector2 const & point, std::vector<Vector2> const & nodesList) const
{
	//https://mathworld.wolfram.com/TriangleInterior.html

	Vector2 v1 = Node(1, nodesList) - Node(0, nodesList);
	Vector2 v2 = Node(2, nodesList) - Node(0, nodesList);

	double a = (Determinant(point, v2) - Determinant(Node(0, nodesList), v2)) / Determinant(v1, v2);
	double b = -1.0F * (Determinant(point, v1) - Determinant(Node(0, nodesList), v1)) / Determinant(v1, v2);

	if (a > 0.0F && b > 0.0F && (a + b) < 1.0F)
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

	for (int i = 0; i < 3; i++)
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

bool Triangle::IsInsideCircumcircle(int vertexID, std::vector<Vector2> const & nodesList) const
{
	//https://stackoverflow.com/questions/39984709/how-can-i-check-wether-a-point-is-inside-the-circumcircle-of-3-points

	Vector2 pos = nodesList[vertexID];

	Matrix_f32 matrix(3, 3);
	Vector2 a = Node(0, nodesList);
	Vector2 b = Node(1, nodesList);
	Vector2 c = Node(2, nodesList);

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

void Triangle::UpdateGeometry(int const * vertixIDs, std::vector<Vector2> const & _nodesList)
{
	UpdateGeometry(vertixIDs[0], vertixIDs[1], vertixIDs[2], _nodesList);
}

void Triangle::UpdateGeometry(int const & vertexID1, int const & vertexID2, int const & vertexID3, std::vector<Vector2> const & nodesList)
{
	vertIDs[0] = vertexID1;
	vertIDs[1] = vertexID2;
	vertIDs[2] = vertexID3;

	area = ComputeArea(Node(0, nodesList), Node(1, nodesList), Node(2, nodesList));
	if (area < 0.0f) //points were passed clockwise, flip than and absolute the area.
	{
		int tempVertID = vertIDs[1];
		vertIDs[1] = vertIDs[2];
		vertIDs[2] = tempVertID;
		area = abs(area); //no need to compute area again, vertex reording only changes sign.
	}
}

Vector2 Triangle::Node(int internalVertexID, std::vector<Vector2> const & nodesList) const
{
	return nodesList[vertIDs[internalVertexID]];
}

Vector2 Triangle::Centroid(std::vector<Vector2> const & nodesList) const
{
	Vector2 centroid;

	centroid = Node(0, nodesList) + Node(1, nodesList) + Node(2, nodesList);
	centroid = Vector2(centroid.x / 3.0, centroid.y / 3.0);

	return centroid;
}

void Triangle::DebugPrintDetails()
{
	std::cout << "-ID: " << id << ".  \tVerts" << vertIDs[0] << ",\t" << vertIDs[1] << ",\t" << vertIDs[2] << std::endl;
}

double Triangle::Determinant(Vector2 const & a, Vector2 const & b) const
{
	return ((static_cast<double>(a.x) * static_cast<double>(b.y)) - (static_cast<double>(a.y) * static_cast<double>(b.x)));
}

double Triangle::ComputeArea(Vector2 node1, Vector2 node2, Vector2 node3) const
{
	return	0.5F *	( static_cast<double>(node1.x) * (static_cast<double>(node2.y) - static_cast<double>(node3.y))
					+ static_cast<double>(node2.x) * (static_cast<double>(node3.y) - static_cast<double>(node1.y))
					+ static_cast<double>(node3.x) * (static_cast<double>(node1.y) - static_cast<double>(node2.y)));
}