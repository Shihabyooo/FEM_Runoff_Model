#include "MeshGenerator.hpp"
#define DICTIONARY_PTR(x) std::unordered_map<size_t, x> *

bool MeshGen::ValidateParameters(MeshGeneratorParameters const & params)
{
	bool status = true;


	//TODO reimplement this.


	return status;
}

bool GriddedTriangulation(MeshGeneratorParameters const & params, std::unordered_map<size_t, Triangle> & outTriList, std::vector<Vector2D> * outNodes, std::vector<size_t> * outBoundaryNodes)
{
	outTriList.clear();

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

		outTriList.insert({ idCounter, subTri1 });
		outTriList.insert({ idCounter + 1, subTri2 });

		idCounter += 2;
	}

	return true;
}

bool MeshGen::GenerateMesh(MeshGeneratorParameters const & params, std::unordered_map<size_t, Triangle> & outTriangles, ::std::vector<Vector2D> * outNodes, ::std::vector<size_t> * outBoundaryNodes)
{
	return GriddedTriangulation(params, outTriangles, outNodes, outBoundaryNodes);
}
