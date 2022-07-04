#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

#include "Globals.hpp"
#include "Triangle.hpp"
#include "LogManager.hpp"

#define MIN_NODES_TO_TRIANGULATE 3
#define SUPER_TRIANGLE_PADDING 0.001 //as percentage of delta of boundary

bool Triangulate(	std::vector<Vector2D> const & nodesList,
					double superTrianglePadding,
					std::unordered_map<int, Triangle> * outTrianglesList,
					std::vector<int> * outBoundaryNodes,
					Triangle * outSuperTriangles);