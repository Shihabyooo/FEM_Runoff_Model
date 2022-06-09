#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

#include "Globals.hpp"
#include "Triangle.hpp"

//TODO bounding box computation should be offloaded to a more global class, since it's needed for rendering.


#define MIN_NODES_TO_TRIANGULATE 3
#define SUPER_TRIANGLE_PADDING 5.0f

//bool Triangulate(std::vector<Vector3> const &nodesList);
bool Triangulate(std::vector<Vector2> nodesList, std::unordered_map<int, Triangle> * outTrianglesList, std::vector<int> * outBoundaryNodes); //Taking a copy of the nodes so I can modify it with the super tri points.