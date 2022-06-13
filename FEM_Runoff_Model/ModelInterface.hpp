#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <GeoTIFF_Parser.h>

#include "Globals.hpp"
#include "DelauneyTriangulation.hpp"

extern std::unordered_map<int, Triangle> triangles;
extern std::vector<Vector2> nodes;
extern std::vector<int> boundaryNodes;
extern Vector2 nodesSW, nodesNE;


bool TestLoadDEM(std::string const & path);
void TestSimulate(std::string const & nodesPath);