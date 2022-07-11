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

	virtual bool Validate()
	{
		return (area > 0.0); //Consider having a minimum viable area for elements tested in  these classes.
	};

	void DebugPrintDetails()
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
	
	//TODO the value bellow is better stored in a some object in ModelInterface. Leaving it here until I fix the element ID ordering to\
	be sequential starting from 0.
	double elementPrecipitation = 0.0; //in m/hr. The precipitation at last pass.	
protected:
	
	int vertCount = 3;
	size_t * vertIDs = NULL;

	std::vector<Vector2D> const * const nodesList;
	double area = 0.0;
};