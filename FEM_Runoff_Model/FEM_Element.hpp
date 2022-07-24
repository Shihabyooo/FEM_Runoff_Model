#pragma once
#include "Globals.hpp"

class Element //useless on its own.
{
public:
	Element() : nodesList(NULL)
	{
	};

	Element(std::vector<Vector2D> const * const _nodesList) : nodesList(_nodesList)
	{
	}
	
	Element(Element const & elem2) : nodesList(elem2.nodesList)
	{
		vertCount = elem2.vertCount;
		if (vertCount > 0 && elem2.vertIDs != NULL)
		{
			vertIDs = new size_t[vertCount];
			for (int i = 0; i < vertCount; i++)
				vertIDs[i] = elem2.vertIDs[i];
		}

		id = elem2.id;
		area = elem2.area;
	}

	~Element() 
	{
		if (vertIDs != NULL)
			delete[] vertIDs;
	};

	

	Vector2D Node(size_t internalVertID) const
	{
		//Leave bounds check to user...
		//return (internalVertID < vertCount && nodesList != NULL ? (*nodesList)[vertIDs[internalVertID]] : Vector2D());
		return (*nodesList)[vertIDs[internalVertID]];
	};

	double Area() const
	{
		return area;
	}

	size_t VertexID(int internalVertID) const
	{
		//Leave bounds check to user...
		//return (internalVertID < vertCount && nodesList != NULL ? vertIDs[internalVertID] : 0);
		return vertIDs[internalVertID];
	};

	size_t const * VertexIDs() const
	{
		return vertIDs;
	}

	int NodeCount() const
	{
		return vertCount;
	}

	bool ContainsVertex(size_t vertexID) const
	{
		bool containsVert = false;

		for (int i = 0; i < vertCount; i++)
			containsVert = containsVert || (vertIDs[i] == vertexID);

		return containsVert;
	}

	virtual bool ContainsPoint(Vector2D const & point) const
	{
		return false;
	}

	bool ContainsEdge(size_t vertID1, size_t vertID2) const //Note: Does not test whether the verts actually make an edge, rather that the two verts are part of element.
	{
		bool check1 = false, check2 = false;

		for (int i = 0; i < vertCount; i++)
		{
			check1 = check1 || (vertIDs[i] == vertID1);
			check2 = check2 || (vertIDs[i] == vertID2);
		}

		return (check1 && check2);
	}

	virtual bool Validate()
	{
		area = ComputeArea();
		return (area > 0.0); //Consider having a minimum viable area for elements tested in  these classes.
	};

	virtual double ComputeArea() const
	{
		return 0.0;
	}

	Vector2D Centroid() const
	{
		Vector2D centroid;

		for (int i = 0; i < vertCount; i++)
			centroid = centroid + Node(i);

		centroid.x *= 1.0 / static_cast<double>(vertCount);
		centroid.y *= 1.0 / static_cast<double>(vertCount);

		return centroid;
	}

	std::vector<size_t> GetOtherVertices(size_t vertexID) //returns other vertices in counter clockwise order
	{
		std::vector<size_t> otherPts;

		int internalVertID = -1;
		for (int i = 0; i < vertCount; i++)
			if (vertIDs[i] == vertexID)
			{
				internalVertID = i;
				break;
			}

		if (internalVertID < 0) //vert not part of this element.
			return otherPts;

		while (otherPts.size() < vertCount - 1)
		{
			internalVertID = internalVertID >= (vertCount - 1) ? 0 : internalVertID + 1;
			otherPts.push_back(vertIDs[internalVertID]);
		}

		return otherPts;
	}

	
	void DebugPrintDetails() const
	{
		std::cout << "Element id: " << id << " - Area: " << area << " - Verts: ";

		if (vertIDs == NULL)
			std::cout << "[DUMMY ELEMENTS]";
		else
			for (int i = 0; i < vertCount; i++)
				std::cout << vertIDs[i] << (i == vertCount - 1 ? "" : ", ");
		std::cout << std::endl;
	}

public:
	size_t id; //of the element itself	

protected: 

	int vertCount = 0;
	size_t * vertIDs = NULL;

	std::vector<Vector2D> const * const nodesList;
	double area = 0.0;
};