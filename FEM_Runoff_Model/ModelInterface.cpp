#include "ModelInterface.hpp"
#include "FileIO.hpp"
#include "Solvers.hpp"

//std::unordered_map<int, Triangle> triangles;
//std::vector<Vector2> nodes;
std::unordered_map<int, Triangle> triangles;
std::vector<Vector2> nodes;
std::vector<int> boundaryNodes;
Vector2 nodesSW, nodesNE;

//Matrices and Vectors
Matrix_f32 globalC, globalPsiX, globalPsiY;

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

void ConstructGlobalConductanceMatrices()
{
	//Psi-X and Psi-Y (for each element) are 3x3 matrices.
	//[Psi-X_e] = 1/6 * delta T *	|	yj-yk	yk-yi	yi-yk	|
	//								|	yj-yk	yk-yi	yi-yk	|
	//								|	yj-yk	yk-yi	yi-yk	|

	//[Psi-Y_e] = 1/6 * delta T *	|	xk-xj	xi-xk	xk-xi	|
	//								|	xk-xj	xi-xk	xk-xi	|
	//								|	xk-xj	xi-xk	xk-xi	|
	//Where x(i/j/k) and y(i/j/k) are the x, y coord of the i/j/kth node.

	//Create global  Psi-X and Psi-y of size (node x nodes), zero initial value.
	//loop over each triangle in mesh
		//loop over permutations of vertIDs (3 verts per triangle = 9 permutations)
			//Psi-X[permutation] += Psi-X_e[localized permutation], e.g. if triangle is verts 3, 5, 8  (i, j, k)\
			Psi-x[3][3] += Psi-x_e[0][0] = 1/6 * dT * (yj-yk), and\
			Psi-x[5][8] += Psi-x_e[1][2] = 1/6 * dT * (yK-yI), etc
			//Ditto for Psi-Y

}

void ConstructGlobalCapacitanceMatrix()
{
	//Capacitance matrix for each element is 3x3 matrix
	//[C_e] = A/3 * |	1	0	0	|
	//				|	0	1	0	|
	//				|	0	0	1	|
	//Where A is the area of element.

	//Construct global matrix similar to procedure in ConstructGlobalConductanceMatrices()
}

void ConstructGlobalBetaMatrix()
{
	//TODO figure out how to construct this matrix (or vector?)

	//Beta matrix for each element is a 3x1 vector
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.

	////Construct global matrix similar to proceudre in ConstructGlobalConductanceMatrices()
	//Nope. This isn't correct...	
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
