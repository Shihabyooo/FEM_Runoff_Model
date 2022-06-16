#include "ModelInterface.hpp"
#include "FileIO.hpp"

//std::unordered_map<int, Triangle> triangles;
//std::vector<Vector2> nodes;
std::unordered_map<int, Triangle> triangles;
std::vector<Vector2> nodes;
std::vector<int> boundaryNodes;
Vector2 nodesSW, nodesNE;

void TestBoundingBox()
{
	//std::cout << "\n computing bounding box in model interface\n";
	nodesSW = nodes[0];
	nodesNE = nodes[0];

	for (auto it = nodes.begin(); it < nodes.end(); it++)
	{
		Print(*it);
		nodesSW.x = Min(it->x, nodesSW.x);
		nodesSW.y = Min(it->y, nodesSW.y);
		nodesNE.x = Max(it->x, nodesNE.x);
		nodesNE.y = Max(it->y, nodesNE.y);
	}

	LogMan::Log("Loaded nodes with bounds: "
				+ std::to_string(nodesSW.x) + ", " + std::to_string(nodesSW.y) + " and "
				+ std::to_string(nodesNE.x) + ", " + std::to_string(nodesNE.y));
}

void TestSimulate(std::string const & nodesPath)
{
	std::cout << "\n test simulation start\n";
	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	LoadCoordinatePairsCSV(nodesPath, nodes);
	
	TestBoundingBox();
	Triangulate(nodes, &triangles, &boundaryNodes);
}
