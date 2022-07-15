#pragma once
#include "Globals.hpp"
#include "DelauneyTriangulation.hpp"
#include "GridMesh.hpp"

//TODO abstract delauney triangluation and GridMesh here

namespace MeshGen {
	bool ValidateParameters(MeshGeneratorParameters const & params);

	bool GenerateMesh(	MeshGeneratorParameters const & params,
						void * outElements,
						::std::vector<Vector2D> * outNodes,
						::std::vector<size_t> * outBoundaryNodes);
}