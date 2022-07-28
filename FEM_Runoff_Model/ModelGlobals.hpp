#pragma once
#include "FileIO.hpp"
#include "MeshGenerator.hpp"
#include "Solvers.hpp"
#include "LogManager.hpp"

extern std::unordered_map<size_t, Triangle> triangles;
extern std::vector<Vector2D> nodes;
extern std::vector<size_t> boundaryNodes;
extern Vector2D nodesSW, nodesNE;
extern Vector2D shedSW, shedNE;
extern std::vector<Vector2D> shedBoundary;

void ComputeBoundingBox(std::vector<Vector2D> const & points, Vector2D & min, Vector2D & max);
bool IsBoundaryNode(size_t nodeID);
bool CheckParameters(ModelParameters const & params);