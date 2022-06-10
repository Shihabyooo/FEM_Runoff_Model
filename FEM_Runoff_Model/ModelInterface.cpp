#include "ModelInterface.hpp"

//std::unordered_map<int, Triangle> triangles;
//std::vector<Vector2> nodes;
std::unordered_map<int, Triangle> triangles;
std::vector<Vector2> nodes;
std::vector<int> boundaryNodes;
Vector2 nodesSW, nodesNE;

void TestBoundingBox()
{
	/*std::cout << "\n computing bounding box in model interface\n";
	nodesSW = nodes[0];
	nodesNE = nodes[0];*/

	for (auto it = nodes.begin(); it < nodes.end(); it++)
	{
		nodesSW.x = Min(it->x, nodesSW.x);
		nodesSW.y = Min(it->y, nodesSW.y);
		nodesNE.x = Max(it->x, nodesNE.x);
		nodesNE.y = Max(it->y, nodesNE.y);
	}

	/*std::cout << "SW corner: ";
	Print(nodesSW);
	std::cout << "NE corner: ";
	Print(nodesNE);*/
}

void TestSimulate(std::string & nodesPath)
{
	std::cout << "\n test simulation start\n";
	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();

	std::ifstream file;
	file.open(nodesPath);
	if (!file.is_open())
	{
		std::cout << "Failed to open " << nodesPath.c_str() << std::endl;
		return;
	}

	char cBuffer = ' ';
	std::string sBuffer = "";
	float cachedFloat = 0;

	while (!file.eof())
	{	
		file.read(&cBuffer, sizeof(cBuffer));
		
		if (cBuffer == ',')
		{
			cachedFloat = atof(sBuffer.c_str());
			sBuffer = "";
		}
		else if (cBuffer == '\n' || file.eof())
		{
			nodes.push_back(Vector2(cachedFloat, atof(sBuffer.c_str())));
			sBuffer = "";
		}
		else
		{
			sBuffer += cBuffer;
		}
	}

	//std::cout << "\n======================================\nLoaded nodes\n";
	//for (auto it = nodes.begin(); it != nodes.end(); ++it)
	//	std::cout << it->x << ", " << it->y << std::endl;
	

	TestBoundingBox();
	Triangulate(nodes, &triangles, &boundaryNodes);

	/*for (auto it = triangles.begin(); it != triangles.end(); ++it)
		std::cout << "tri: " << it->second.id << " - verts:" << it->second.vertIDs[0] << ", " << it->second.vertIDs[1] << ", " << it->second.vertIDs[2] << std::endl;*/

	//std::cout << "Count: " << nodes.size() << std::endl;
}
