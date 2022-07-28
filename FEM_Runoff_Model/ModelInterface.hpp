#pragma once
#include "Globals.hpp"
#include "MeshGenerator.hpp"

bool LoadWatershedBoundary(std::string const & boundaryPath);
bool GenerateMesh(MeshGeneratorParameters & params);
Triangle const * GetElementContainingPoint(Vector2D const & pos);
bool UpdateNode(size_t id, Vector2D const & newPos);

bool LoadTimeSeries(std::string const & path, TimeSeries & ts);
bool Simulate(ModelParameters const & params);

//Getters (mostly for GUI)
std::vector<Vector2D> const & GetNodes();
std::unordered_map<size_t, Triangle> const & GetTriangles();
std::vector<size_t> const & GetBoundaryNodes();
std::pair<Vector2D const &, Vector2D const &> GetNodesBoundingBox();
std::pair<Vector2D const &, Vector2D const &> GetWatershedBoundingBox();
std::vector<Vector2D> const & GetWatershedBoundary();