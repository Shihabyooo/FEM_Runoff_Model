#include "Triangle.hpp"

Triangle::Triangle(int _id, int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2> const & nodesList)
{
	id = _id;
	//nodesList = _nodesList;
	UpdateGeometry(vertexID1, vertexID2, vertexID3, nodesList);
}

Triangle::~Triangle()
{
}

void Triangle::Subdivide(int centroidID, int baseTriangleID, std::vector<Vector2> const & nodesList, std::vector<Triangle>& outNewTriangles)
{
	//Subdivision stage

	outNewTriangles.push_back(Triangle(baseTriangleID, centroidID, vertIDs[0], vertIDs[1], nodesList));
	outNewTriangles.push_back(Triangle(baseTriangleID + 1, centroidID, vertIDs[0], vertIDs[2], nodesList));
	outNewTriangles.push_back(Triangle(baseTriangleID + 2, centroidID, vertIDs[1], vertIDs[2], nodesList));

	//Optimization stage. If the current centroid is within the cirumscribed circle (circumcircle) of the neighbouring triangle of each of our new 
	//triangles, flip the joint edge between them so it connects the non connected vertices (the centroid and the third vert of neighbour) instead.
}

bool Triangle::ContainsPoint(Vector2 const & point, std::vector<Vector2> const & nodesList)
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
	return (vertIDs[0] == vertexID) || (vertIDs[1] == vertexID) || (vertIDs[2] == vertexID);
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

bool Triangle::IsExternalTriangle(int * externalVertices, int * outMeshEdgeVerts, int * outMeshEdgeContribCount)
{
	*outMeshEdgeContribCount = 0;
	
	for (int i = 0; i < 3; i++)
	{
		bool isExternalVert = false;
		for (int j = 0; j < 3; j++)
		{
			if (externalVertices[j] == vertIDs[i])
			{
				isExternalVert = true;
				break;
			}
		}
		if (!isExternalVert)
		{
			if (*outMeshEdgeContribCount >= 2)
				return false;

			outMeshEdgeVerts[*outMeshEdgeContribCount] = vertIDs[i];
			*outMeshEdgeContribCount++;
		}
	}

	return (*outMeshEdgeContribCount >= 2);
}

bool Triangle::IsNeighbour(Triangle const & triangle, int outsideVertex, int ** outSharedVertsIDs)
{
	std::vector<int> edgeVerts;

	for (int i = 0; i < 3; i++)
	{
		if (vertIDs[i] != outsideVertex)
			edgeVerts.push_back(vertIDs[i]);
	}

	if (triangle.ContainsEdge(edgeVerts[0], edgeVerts[1]))
	{
		*outSharedVertsIDs = new int[2] {edgeVerts[0], edgeVerts[1]};
		return true;
	}
	else
	{
		*outSharedVertsIDs = NULL;
		return false;
	}
}

bool Triangle::IsInsideCircumcircle(int vertexID, std::vector<Vector2> const & nodesList)
{
	//https://stackoverflow.com/questions/39984709/how-can-i-check-wether-a-point-is-inside-the-circumcircle-of-3-points

	Vector2 pos = nodesList[vertexID];

	Matrix_f32 matrix(4, 4);
	matrix[0][0] = Node(0, nodesList).x;
	matrix[0][1] = Node(0, nodesList).y;
	matrix[0][2] = Node(0, nodesList).x * Node(0, nodesList).x + Node(0, nodesList).y * Node(0, nodesList).y;
	matrix[0][3] = 1.0f;

	matrix[1][0] = Node(1, nodesList).x;
	matrix[1][1] = Node(1, nodesList).y;
	matrix[1][2] = Node(1, nodesList).x * Node(1, nodesList).x + Node(1, nodesList).y * Node(1, nodesList).y;
	matrix[1][3] = 1.0f;

	matrix[2][0] = Node(2, nodesList).x;
	matrix[2][1] = Node(2, nodesList).y;
	matrix[2][2] = Node(2, nodesList).x * Node(2, nodesList).x + Node(2, nodesList).y * Node(2, nodesList).y;
	matrix[2][3] = 1.0f;

	matrix[3][0] = pos.x;
	matrix[3][1] = pos.y;
	matrix[3][2] = pos.x * pos.x + pos.y * pos.y;
	matrix[3][3] = 1.0f;

	return (matrix.Determinant() > 0.0f);
}

void Triangle::UpdateGeometry(int vertexID1, int vertexID2, int vertexID3, std::vector<Vector2> const & nodesList)
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

double Triangle::Determinant(Vector2 const & a, Vector2 const & b)
{
	return (a.x * b.y) - (a.y * b.x);
}

double Triangle::ComputeArea(Vector2 node1, Vector2 node2, Vector2 node3)
{
	return	0.5F * node1.x * (node2.y - node3.y) + node2.x * (node3.y - node1.y) + node3.x * (node1.y - node2.y);
}