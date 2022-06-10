#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "Globals.hpp"
#include "DelauneyTriangulation.hpp"

extern std::unordered_map<int, Triangle> triangles;
extern std::vector<Vector2> nodes;
extern std::vector<int> boundaryNodes;
extern Vector2 nodesSW, nodesNE;

void TestSimulate(std::string & nodesPath);