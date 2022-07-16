#include "MeshGenerator.hpp"
#define DICTIONARY_PTR(x) std::unordered_map<size_t, x> *

bool MeshGen::ValidateParameters(MeshGeneratorParameters const & params)
{
	bool status = true;

	if (params.meshType == ElementType::undefined)
	{
		LogMan::Log("ERROR! Mesh type must be defined!", LOG_WARN);
		status = false;
	}
	else if (params.meshType == ElementType::triangle)
	{
		/*if (params.useCustomNodes && params.inNodesList == NULL)
		{
			LogMan::Log("ERROR! Input nodes must be set when using custom nodes for triangulation!", LOG_WARN);
			status = false;
		}*/
		if (params.superTrianglePadding <= 0.0)
		{
			LogMan::Log("ERROR! Super triangle padding must be greater than 0.0!", LOG_WARN);
			status = false;
		}
		
		if (!params.useCustomNodes && params.resolution < 2)
		{
			LogMan::Log("ERROR! Resolution cannot be less than 2!", LOG_WARN);
			status = false;
		}
	}
	else if (params.meshType == ElementType::rectangle)
	{
		if (params.resolution < 2)
		{
			LogMan::Log("ERROR! Resolution cannot be less than 2!", LOG_WARN);
			status = false;
		}
		if (params.internalPadding <= 0.0)
		{
			LogMan::Log("ERROR! Internal grid padding must be greater than 0.0!", LOG_WARN);
			status = false;
		}
		if (params.rayCastPadding <= 0.0)
		{
			LogMan::Log("ERROR! Raycast padding must be greater than 0.0!", LOG_WARN);
			status = false;
		}
	}

	return status;
}

bool GriddedTriangulation(MeshGeneratorParameters const & params, std::unordered_map<size_t, Triangle> * outTriList, std::vector<Vector2D> * outNodes, std::vector<size_t> * outBoundaryNodes)
{
	//generate a temporary rectangular grid, then subdivide it into triangles.
	std::unordered_map<size_t, Rectangle> tempRectElements;
	bool status = GenerateGrid(*params.boundary, outNodes, &tempRectElements, outBoundaryNodes, params.resolution, params.internalPadding, params.rayCastPadding);

	if (!status)
		return false;

	size_t idCounter = 0;
	for (auto it = tempRectElements.begin(); it != tempRectElements.end(); ++it)
	{
		size_t const * vertices = it->second.VertexIDs();

		Triangle subTri1(idCounter, vertices[0], vertices[1], vertices[3], outNodes);
		Triangle subTri2(idCounter + 1, vertices[1], vertices[2], vertices[3], outNodes);

		outTriList->insert({ idCounter, subTri1 });
		outTriList->insert({ idCounter + 1, subTri2 });

		idCounter += 2;
	}

	return true;
}

bool MeshGen::GenerateMesh(MeshGeneratorParameters const & params, void * outElements, ::std::vector<Vector2D> * outNodes, ::std::vector<size_t> * outBoundaryNodes)
{
	switch (params.meshType)
	{
	case ElementType::rectangle:
	{
		std::unordered_map<size_t, Rectangle> * rectOutput = static_cast<std::unordered_map<size_t, Rectangle> *>(outElements);
		return GenerateGrid(*params.boundary, outNodes, rectOutput, outBoundaryNodes, params.resolution, params.internalPadding, params.rayCastPadding);
	}
	case ElementType::triangle:
	{
		std::unordered_map<size_t, Triangle> * triOutput = static_cast<std::unordered_map<size_t, Triangle> *>(outElements);
		if (params.useCustomNodes)
			return Triangulate(*(params.inNodesList), params.superTrianglePadding, triOutput, outBoundaryNodes, params.outSuperTriangleNodes);
		else
			return GriddedTriangulation(params, triOutput, outNodes, outBoundaryNodes);
	}
	default:
		LogMan::Log("ERROR! Internal error! (caught invalid ElementType in GenerateMesh())", LOG_ERROR);
		return false;
	}
}
