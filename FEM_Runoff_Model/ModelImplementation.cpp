#pragma once
#include "ModelImplementation.hpp"


//TODO add a cleanup method to clear the allocated memory (superTriangles, rasters) when program closes.
ElementType activeMeshType = ElementType::undefined;
std::unordered_map<size_t, Rectangle> rectangles;
std::unordered_map<size_t, Triangle> triangles;
Vector2D superTriangles[6];
std::vector<Vector2D> nodes;
std::vector<size_t> boundaryNodes;
Vector2D nodesSW, nodesNE;
Vector2D shedSW, shedNE;
std::vector<Vector2D> shedBoundary;

Vector_f64 heads;

//Precipitation returned as meters per hour.
//current impl doesn't need triangle, but later it would.
double GetCurrentPrecipitation(double time, ModelParameters const & params, Element const * element)
{
	if (params.variablePrecipitation)
	{
		LogMan::Log("Warning! Variable precipitation is not yet implemented!", LOG_WARN);
		return 0.0;
	}
	else
	{
		return params.unitTimeSeries.SampleRate(time) / 1000.0;// , params.timeStep, params.precipitationTemporalInterpolationType);
	}
}

Vector_f64 ComputePreciptationVector(double time, ModelParameters const & params)
{
	//Preciptation Vector is the last term of the RHS of the formulation. i.e. dT * {Beta} * ((1 - omega) * Pe_t + omega * Pe_t+dt)

	//Beta matrix for each element is a 3x1 vector
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.

	//zero out outVector before doing anything.
	Vector_f64 result(nodes.size());

	size_t elementCount = 0;
	int nodesPerElement = 0;

	switch (activeMeshType)
	{
	case ElementType::rectangle:
		elementCount = rectangles.size();
		nodesPerElement = 4;
		break;
	case ElementType::triangle:
		elementCount = triangles.size();
		nodesPerElement = 3;
		break;
	}

	for (size_t i = 0; i < elementCount; i++)
	{
		Element * element = NULL;

		switch (activeMeshType)
		{
		case ElementType::rectangle:
			element = new Element(rectangles[i]);
			break;
		case ElementType::triangle:
			element = new Element(triangles[i]);
			break;
		}

		double newPrecipitation = GetCurrentPrecipitation(time + params.timeStep, params, element);
		double elementContrib = (element->Area() / 3.0) * ((1.0 - params.femOmega) * element->elementPrecipitation + params.femOmega * newPrecipitation);
		elementContrib *= params.timeStep;
		element->elementPrecipitation = newPrecipitation;

		//test
		/*if (!element->ContainsVertex(1) && !element->ContainsVertex(2) && !element->ContainsVertex(3))
			elementContrib = 0.0;*/
			//end test

		for (int intNodeID = 0; intNodeID < element->NodeCount(); intNodeID++)
			result[element->VertexID(intNodeID)] += elementContrib;

		delete element;
	}

	return result;
}

//TODO have to vector_f64 to store computed nodal velocities, and a function to compute the newheads nodal velocities every \
iteration. newHead velos are then moved to oldHead velos after internal iterations are done.
std::pair<double, double> ComputeVelocityComponents(size_t nodeID, Vector_f64 waterElevation) //waterElevation = nodeElevation + head
{
	//if (IsBoundaryNode(nodeID)) //test
	//	return std::pair<double, double>(0.0, 0.0); //test

	double head = waterElevation[nodeID] - nodeElevation[nodeID];

	double velocity = 3600.0 * sqrt(nodeSlope[nodeID]) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];

	double angle = FDR2Angle(lround(nodeFDR[nodeID]));

	double u = velocity * cos(angle);
	double v = velocity * sin(angle);

	/*double u = Sign(nodeSlopeX[nodeID]) * 3600.0 * sqrt(abs(nodeSlopeX[nodeID])) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];
	double v = Sign(nodeSlopeY[nodeID]) *3600.0 * sqrt(abs(nodeSlopeY[nodeID])) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];*/

	return std::pair<double, double>(u, v);
}

Matrix_f64 ElementKMatrix_Tri(Triangle const & tri, std::pair<double, double> uv[3])
{
	//						|	ui (yj - yk) + vi (xk - xj)		ui (yk - yi) + vi (xi - xk)		ui (yi - yj) + vi (xj - xi)	|
	//[K_element] = 1/6	*	|	uj (yj - yk) + vj (xk - xj)		uj (yk - yi) + vj (xi - xk)		uj (yi - yj) + vj (xj - xi)	|
	//						|	uk (yj - yk) + vk (xk - xj)		uk (yk - yi) + vk (xi - xk)		uk (yi - yj) + vk (xj - xi)	|

	//The 1/6 scalar is mulitplied in the calling function.

	Vector2D const & i = tri.Node(0);
	Vector2D const & j = tri.Node(1);
	Vector2D const & k = tri.Node(2);

	Matrix_f64 kMat(3, 3);

	//uv[0].second = uv[1].second = uv[2].second = (uv[0].second + uv[1].second + uv[2].second) / 3.0;//test
	//uv[0].first = uv[1].first = uv[2].first = (uv[0].first + uv[2].first + uv[2].first) / 3.0;//test

	kMat[0][0] = uv[0].first * (j.y - k.y) + uv[0].second * (k.x - j.x);
	kMat[0][1] = uv[0].first * (k.y - i.y) + uv[0].second * (i.x - k.x);
	kMat[0][2] = uv[0].first * (i.y - j.y) + uv[0].second * (j.x - i.x);

	kMat[1][0] = uv[1].first * (j.y - k.y) + uv[1].second * (k.x - j.x);
	kMat[1][1] = uv[1].first * (k.y - i.y) + uv[1].second * (i.x - k.x);
	kMat[1][2] = uv[1].first * (i.y - j.y) + uv[1].second * (j.x - i.x);

	kMat[2][0] = uv[2].first * (j.y - k.y) + uv[2].second * (k.x - j.x);
	kMat[2][1] = uv[2].first * (k.y - i.y) + uv[2].second * (i.x - k.x);
	kMat[2][2] = uv[2].first * (i.y - j.y) + uv[2].second * (j.x - i.x);

	return kMat;
}

double k_e_x_fixed[4][4]{ {-2, 2, 1, -1},
							{-2, 2, 1, -1},
							{-1, 1, 2, -2},
							{-1, 1, 2, -2} };

double k_e_y_fixed[4][4]{ {-2, -1, 1, 2},
							{-1, -2, 2, 1},
							{-1, -2, 2, 1},
							{-2, -1, 1, 2} };

Matrix_f64 ElementKMatrix_Rect(Rectangle const & rect, std::pair<double, double> uv[4])
{
	//[K_element] =	[K_e_x] + [K_e_y]
	//
	//								|	ui	0	0	0	|		|	-2	2	1	-1	|
	//[K_e_x]	=	height/12.0 *	|	0	uj	0	0	|	*	|	-2	2	1	-1	|
	//								|	0	0	uk	0	|		|	-1	1	2	-2	|
	//								|	0	0	0	ul	|		|	-1	1	2	-2	|
	//
	//								|	vi	0	0	0	|		|	-2	-1	1	2	|
	//[K_e_y]	=	width/12.0 *	|	0	vj	0	0	|	*	|	-1	-2	2	1	|
	//								|	0	0	vk	0	|		|	-1	-2	2	1	|
	//								|	0	0	0	vl	|		|	-2	-1	1	2	|
	//

	Matrix_f64 kMat_x(4, 4);
	Matrix_f64 kMat_y(4, 4);

	double h12 = rect.Height() / 12.0;
	double w12 = rect.Width() / 12.0;
	double multX, multY;

	for (int i = 0; i < 4; i++)
	{
		multX = uv[i].first * h12;
		multY = uv[i].second * w12;

		for (int j = 0; j < 4; j++)
		{
			kMat_x[i][j] = multX * k_e_x_fixed[i][j];
			kMat_y[i][j] = multY * k_e_y_fixed[i][j];
		}
	}

	//std::cout << "\n//////////////////////////////////////////\n";
	//std::cout << "element: "; 
	//rect.DebugPrintDetails();
	//(kMat_x + kMat_y).DisplayOnCLI(2);//test

	return kMat_x + kMat_y;
}

Matrix_f64 ElementCMatrixLumped(int size, double areaDivision)
{
	//Lumped Capcitance Matrix for both element type defined as a diagonal matrix with diag value equal to element area / number of nodes\
	i.e. "areaDivision."

	Matrix_f64 cMat(size, size);

	for (int i = 0; i < size; i++)
		cMat[i][i] = areaDivision;

	return cMat;
}

Matrix_f64 ElementCMatrix_Tri(ModelParameters const & params, Triangle const & tri)
{
	if (params.useLumpedForm)
		return ElementCMatrixLumped(3, tri.Area() / 3.0);

	//else return Constistance capcitance matrix

	Matrix_f64 cMat(3, 3);

	double mult = tri.Area() / 12.0;

	cMat[0][0] = cMat[1][1] = cMat[2][2] = mult * 2.0;
	cMat[0][1] = cMat[0][2] = cMat[1][0] = cMat[1][2] = cMat[2][0] = cMat[2][1] = mult;

	return cMat;
}

double c_e_fixed[4][4]{ {4, 2, 1, 2},
						{2, 4, 2, 1},
						{1, 2, 4, 2},
						{2, 1, 2, 4} };

Matrix_f64 ElementCMatrix_Rect(ModelParameters const & params, Rectangle const & rect)
{
	if (params.useLumpedForm)
		return ElementCMatrixLumped(4, rect.Area() / 4.0);

	//else return Constistance capcitance matrix
	Matrix_f64 cMat(4, 4);
	double mult = rect.Area() / 36.0;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			cMat[i][j] = c_e_fixed[i][j] * mult;

	return cMat;
}

Matrix_f64 ComputeGlobalCoefficientsMatrix(ModelParameters const & params, Vector_f64 const & newHeads)
{
	//[C] + w*dt*[K]
	Matrix_f64 coefMat(nodes.size(), nodes.size());
	Matrix_f64 kMat, cMat;

	if (activeMeshType == ElementType::triangle)
	{
		for (auto it = triangles.begin(); it != triangles.end(); ++it)
		{
			std::pair<double, double> uv[3] = { ComputeVelocityComponents(it->second.VertexID(0), newHeads),
												ComputeVelocityComponents(it->second.VertexID(1), newHeads),
												ComputeVelocityComponents(it->second.VertexID(2), newHeads) };

			kMat = ElementKMatrix_Tri(it->second, uv) * (params.timeStep * params.femOmega / 6.0);
			cMat = ElementCMatrix_Tri(params, it->second);

			for (int row = 0; row < 3; row++)
				for (int column = 0; column < 3; column++)
					coefMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] + kMat[row][column]);
		}
	}
	else
	{
		for (auto it = rectangles.begin(); it != rectangles.end(); ++it)
		{
			std::pair<double, double> uv[4] = { ComputeVelocityComponents(it->second.VertexID(0), newHeads),
												ComputeVelocityComponents(it->second.VertexID(1), newHeads),
												ComputeVelocityComponents(it->second.VertexID(2), newHeads),
												ComputeVelocityComponents(it->second.VertexID(3), newHeads) };

			kMat = ElementKMatrix_Rect(it->second, uv) * (params.timeStep * params.femOmega / 6.0);
			cMat = ElementCMatrix_Rect(params, it->second);

			cMat[0][0] = cMat[1][1] = cMat[2][2] = cMat[3][3] = (it->second.Area() / 4.0);

			for (int row = 0; row < 4; row++)
				for (int column = 0; column < 4; column++)
					coefMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] + kMat[row][column]);

			//test
			/*if (it->second.ContainsVertex(46))
			{
				std::cout << "Inside CoefMatConst\nElem: ";
				it->second.DebugPrintDetails();
				kMat.DisplayOnCLI();
			}*/
		}
	}

	return coefMat;
}

Matrix_f64 ComputeGlobalConductanceMatrix(ModelParameters const & params)
{
	//[C] - dt*(1-w)*[K]

	Matrix_f64 condMat(nodes.size(), nodes.size());
	Matrix_f64 kMat;

	if (activeMeshType == ElementType::triangle)
	{
		Matrix_f64 cMat(3, 3);

		for (auto it = triangles.begin(); it != triangles.end(); ++it)
		{
			std::pair<double, double> uv[3] = { ComputeVelocityComponents(it->second.VertexID(0), heads),
												ComputeVelocityComponents(it->second.VertexID(1), heads),
												ComputeVelocityComponents(it->second.VertexID(2), heads) };

			kMat = ElementKMatrix_Tri(it->second, uv) * (params.timeStep * (1.0 - params.femOmega) / 6.0);

			cMat[0][0] = cMat[1][1] = cMat[2][2] = (it->second.Area() / 3.0);

			for (int row = 0; row < 3; row++)
				for (int column = 0; column < 3; column++)
					condMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] - kMat[row][column]);
		}
	}
	else
	{
		Matrix_f64 cMat(4, 4);

		for (auto it = rectangles.begin(); it != rectangles.end(); ++it)
		{
			std::pair<double, double> uv[4] = { ComputeVelocityComponents(it->second.VertexID(0), heads),
												ComputeVelocityComponents(it->second.VertexID(1), heads),
												ComputeVelocityComponents(it->second.VertexID(2), heads),
												ComputeVelocityComponents(it->second.VertexID(3), heads) };

			kMat = ElementKMatrix_Rect(it->second, uv) * (params.timeStep * (1.0 - params.femOmega) / 6.0);

			cMat[0][0] = cMat[1][1] = cMat[2][2] = cMat[3][3] = (it->second.Area() / 4.0);

			for (int row = 0; row < 4; row++)
				for (int column = 0; column < 4; column++)
					condMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] - kMat[row][column]);

			//test
			/*if (it->second.ContainsVertex(46))
			{
				std::cout << "Inside CondMatConst\nElem: ";
				it->second.DebugPrintDetails();
				kMat.DisplayOnCLI();
			}*/

			//test
			/*if (!kMat.Determinant() < 0.00001)
			{
				std::cout << "ping! at elem:";
				it->second.DebugPrintDetails();
			}*/
		}
	}

	return condMat;
}

#pragma region Test

Vector_f64 _new_h, _precipComp, _RHS;
std::vector<std::pair<double, double>> qTS;
void TestShowValues(ModelParameters const & params)
{
	std::cout << "\n  id     |  head | newH |||  |||   u0   |    v0   |   u1    |   v1    |  precip  |  RHS\n";
	for (size_t i = 0; i < nodes.size(); i++)
	{
		auto uv0 = ComputeVelocityComponents(i, heads);
		auto uv1 = ComputeVelocityComponents(i, _new_h);

		std::cout << std::setw(4) << i << " " << (IsBoundaryNode(i) ? "[B]" : "   ") << " | " <<
			std::fixed << std::setprecision(2) << std::setw(5) << heads[i] << " | " << std::setw(5) << _new_h[i] <<
			"|||" << (heads[i] < _new_h[i] ? "UP" : (heads[i] > _new_h[i] ? "DN" : "--")) << "|||" <<
			std::fixed << std::setprecision(2) << std::setw(7) << uv0.first << " | " << std::setw(7) << uv0.second << " | " <<
			std::fixed << std::setprecision(2) << std::setw(7) << uv1.first << " | " << std::setw(7) << uv1.second << " | " <<
			std::fixed << std::setprecision(2) << std::setw(7) << _precipComp[i] << " | " <<
			std::fixed << std::setprecision(2) << std::setw(7) << _RHS[i] << std::endl;
	}


	size_t targetRowToShow = 0;
	std::cout << "\nConductance Matrix:\n";
	Matrix_f64 _globalConduc = ComputeGlobalConductanceMatrix(params);
	//(globalC -  ComputeGlobalConductanceMatrix(params) ).DisplayOnCLI(0);
	//ComputeGlobalConductanceMatrix(params).DisplayOnCLI(0);
	for (size_t i = 0; i < _globalConduc.Columns(); i++)
		std::cout << std::fixed << std::setprecision(0) << std::setw(0) << _globalConduc[targetRowToShow][i] << " ";
	std::cout << "\n";


	std::cout << "\nCoefficients Matrix:\n";
	Matrix_f64 _globalCoef = ComputeGlobalCoefficientsMatrix(params, _new_h);
	//(ComputeGlobalCoefficientsMatrix(params, _new_h) - globalC).DisplayOnCLI(0);
	//ComputeGlobalCoefficientsMatrix(params, _new_h).DisplayOnCLI(0);
	for (size_t i = 0; i < _globalCoef.Columns(); i++)
		std::cout << std::fixed << std::setprecision(0) << std::setw(0) << _globalCoef[targetRowToShow][i] << " ";
	std::cout << "\n";


	std::cout << std::endl;
}
#pragma endregion

void ComputeRHSVector(double time, ModelParameters const & params, Vector_f64 & outRHS)
{
	//[GlobalConductanceMat] * {h_0} + precipComponent

	outRHS = ComputeGlobalConductanceMatrix(params) * heads + ComputePreciptationVector(time, params);

	//test
	_RHS = outRHS;
	_precipComp = ComputePreciptationVector(time, params);

	//TestShowValues(params);
}

//[A]{x} = {b}
void AdjustForBoundaryConditions(Matrix_f64 & aMat, Vector_f64 & xVec, Vector_f64 & bVec)
{
	////Using https://finite-element.github.io/7_boundary_conditions.html
	//for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	//{
	//	bVec[*it] = nodeElevation[*it];

	//	for (size_t i = 0; i < aMat.Columns(); i++)
	//		aMat[*it][i] = 0.0;

	//	aMat[*it][*it] = 1.0;
	//}

	//Using Istok's method
	size_t reducedSystemSize = nodes.size() - boundaryNodes.size();
	Matrix_f64 adjustedAMat(reducedSystemSize, reducedSystemSize);
	Vector_f64 adjustedXVec(reducedSystemSize);
	Vector_f64 adjustedBVec(reducedSystemSize);

	size_t rowCounter = 0;
	size_t columnCounter = 0;
	for (size_t row = 0; row < aMat.Rows(); row++)
	{
		if (IsBoundaryNode(row))
			continue;

		adjustedBVec[rowCounter] = bVec[row];
		//adjustedXVec[rowCounter] = xVec[row]; //pointless. xVec is not filled with data yet.

		for (size_t column = 0; column < aMat.Columns(); column++)
		{
			if (IsBoundaryNode(column))
				continue;

			adjustedAMat[rowCounter][columnCounter] = aMat[row][column];
			columnCounter++;
		}
		columnCounter = 0;
		rowCounter++;
	}

	aMat = adjustedAMat;
	bVec = adjustedBVec;
	xVec = adjustedXVec;
}

void ReIntroduceBoundaryNodes(Vector_f64 & vec)
{
	/*if (vec.Rows() == nodes.size())
		return;*/

	Vector_f64 adjustedVec(nodes.size());

	size_t rowCounter = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (IsBoundaryNode(i))
		{
			adjustedVec[i] = 0.0;
		}
		else
		{
			adjustedVec[i] = vec[rowCounter];
			rowCounter++;
		}
	}

	vec = adjustedVec;
}

bool RunSimulation(ModelParameters const & params)
{
	LogMan::Log("Starting a simulation run");

#pragma region Init
	if (!CheckParameters(params))
		return false;

	UnloadAllRasters();

	if (!LoadInputRasters(params))
		return false;

	//Special consideration. Since the boundary node listing includes our exit node, we have to manually remove it.
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == params.outletNode)
		{
			boundaryNodes.erase(it);
			break;
		}

	LogMan::Log("Caching hydraulic variables.");

	if (!CacheManningCoefficients(params))
		return false;

	if (!CacheSlopes(params))
		return false;

	if (!CacheElevations())
		return false;

#pragma endregion

#pragma region Test
	/*std::cout << "\n===================================================\n";
	std::cout << "Global Capacitance";
	std::cout << "\n===================================================\n";
	globalC.DisplayOnCLI(0);*/

	std::cout << "\n===================================================\n";
	std::cout << "Time series";
	std::cout << "\n===================================================\n";
	for (int i = 0; i < params.unitTimeSeries.size; i++)
		std::cout << params.unitTimeSeries.series[i].first << " - " << params.unitTimeSeries.series[i].second << std::endl;

	std::cout << "\n===================================================\n";
	std::cout << "Nodes";
	std::cout << "\n===================================================\n";
	size_t nodeIDCounter = 0;
	for (auto it = nodes.begin(); it != nodes.end(); ++it)
	{
		std::cout << nodeIDCounter << "\t - \t" << std::fixed;
		Print(*it);
		nodeIDCounter++;
	}

	std::cout << "\n===================================================\n";
	std::cout << "Boundery Nodes";
	std::cout << "\n===================================================\n";
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		std::cout << *it << std::endl;

	/*std::cout << "\n===================================================\n";
	std::cout << "Elevations, Slopes, Manning roughness coef";
	std::cout << "\n===================================================\n";*/
	/*std::cout << "node | Elev  |   n  |  Sx |  Sy |  FDR  \n";
	for (size_t i = 0; i < nodes.size(); i++)
		std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " << nodeElevation[i] << " | " <<
		nodeManning[i] << " | " << nodeSlopeX[i] << " | " << nodeSlopeY[i] << " | "  << nodeFDR[i] << std::endl;*/
		/*std::cout << "node | Elev  |   n  |  Slope  |  FDR  \n";
		for (size_t i = 0; i < nodes.size(); i++)
			std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " << nodeElevation[i] << " | "  <<
						nodeManning[i] << " | " << nodeSlope[i] << " | " << nodeFDR[i]<< std::endl;*/

						//return false;
#pragma endregion

	//test
	//nodeFDR[9] = 4;
	//nodeFDR[15] = 4;
	//nodeFDR[20] = 3;
	//nodeFDR[21] = 3;
	//nodeFDR[22] = 4;
	//nodeFDR[32] = 4;
	//nodeFDR[55] = 3;
	//nodeFDR[56] = 4;
	//nodeFDR[63] = 5;

	//endtest



	nodeElevation = Vector_f64(nodes.size()); //test
	heads = nodeElevation;

	//heads[7] += 2.0;//test
	//heads[57] += 2.0;//test
	//heads[248] += 2.0;//test

	//Loop from start time to end time
	double time = params.startTime;
	LogMan::Log("Starting simulation loop");

	while (time <= params.endTime)
	{
		std::cout << "\n\n===================================================\n";
		LogMan::Log("At T= " + std::to_string(time));
		std::cout << "===================================================\n";

		//Vector_f64 newHeads = heads + Vector_f64(nodes.size(), 0.001);
		Vector_f64 newHeads(nodes.size());
		for (size_t i = 0; i < nodes.size(); i++)
			if (!IsBoundaryNode(i))
				newHeads[i] = heads[i] + 0.001;

		//internal loop
		for (size_t i = 0; i <= params.maxInternalIterations; i++)
		{
			_new_h = newHeads; //test

			Vector_f64 RHS;
			ComputeRHSVector(time, params, RHS);

			Vector_f64 residuals; //needed by solver
			Vector_f64 fixedNewH;

			Matrix_f64 coeffMat = ComputeGlobalCoefficientsMatrix(params, newHeads);

			AdjustForBoundaryConditions(coeffMat, fixedNewH, RHS);

			if (!Solve(coeffMat, RHS, fixedNewH, residuals, params))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}

			//force computed newHeads to be positive (not sure about this)
			for (size_t j = 0; j < fixedNewH.Rows(); j++)
				fixedNewH[j] = Max(fixedNewH[j], nodeElevation[j]);

			ReIntroduceBoundaryNodes(fixedNewH);

			std::cout << "Solver residual: " << std::fixed << std::setprecision(10) << residuals.Magnitude() << std::endl; //test
			std::cout << "current Internal Residual: " << std::fixed << std::setprecision(10) << (newHeads - fixedNewH).Magnitude() << std::endl; //test
			if ((newHeads - fixedNewH).Magnitude() <= params.internalResidualTreshold)
			{
				newHeads = fixedNewH;
				break;
			}
			else if (i >= params.maxInternalIterations)
				LogMan::Log("Reached maxInternalIterations without reaching appropriate h", LOG_WARN);

			newHeads = fixedNewH;
		}
		//TestShowValues(params);

		std::cout << "\n------------------------------------------------------\n";
		std::cout << "heads result at time: " << std::setprecision(4) << time << " --> " << std::setprecision(4) << time + params.timeStep << std::endl;
		double headSum = 0.0, newHeadSum = 0.0;
		for (size_t i = 0; i < heads.Rows(); i++)
		{
			//std::cout << i << "\t slopes: " << std::setprecision(4) << nodeSlopeX[i] << " - " << nodeSlopeY[i] << "\t - \t" << std::setprecision(4) << heads[i] << "\t-->\t" <<
			std::cout << i << "\t slopes: " << std::setprecision(4) << nodeSlope[i] << " - " << std::setprecision(0) << nodeFDR[i] << "\t - \t" << std::setprecision(4) << heads[i] << "\t-->\t" <<
				std::setprecision(4) << newHeads[i] << "  " <<
				(newHeads[i] > heads[i] ? "UP" : (newHeads[i] == heads[i] ? "--" : "DN")) <<
				std::endl;

			headSum += heads[i];
			newHeadSum += newHeads[i];
		}
		std::cout << "=============== Sum : " << headSum << " --> " << newHeadSum << std::endl;
		std::cout << "\n------------------------------------------------------\n";
		/*double qx = sqrt(nodeSlopeX[params.outletNode]) * pow(heads[params.outletNode], 5.0 / 3.0) / nodeManning[params.outletNode];
		double qy = sqrt(nodeSlopeY[params.outletNode]) * pow(heads[params.outletNode], 5.0 / 3.0) / nodeManning[params.outletNode];
		double q = sqrt(qx * qx + qy * qy);*/
		double q = sqrt(nodeSlope[params.outletNode]) * pow(heads[params.outletNode], 5.0 / 3.0) / nodeManning[params.outletNode];
		double flowWidth;
		if (activeMeshType == ElementType::triangle)
			flowWidth = 2.0 *  abs((triangles[0].Centroid() - triangles[0].Node(0)).x); //test. 
		else
			flowWidth = rectangles[0].Width(); //test. 
		std::cout << "h: " << std::setprecision(7) << heads[params.outletNode] << "\tQ: " << std::setprecision(3) << (q * flowWidth) << std::endl;
		qTS.push_back(std::pair<double, double>(time, q * flowWidth));
		std::cout << "\n------------------------------------------------------\n";

		heads = newHeads;
		time += params.timeStep;

		//test
		std::cout << "Enter to proceed to next step\n";
		std::cin.sync();
		std::cin.get();
	}

	//test
	std::cout << " time \t Q (cms)\n";
	for (auto it = qTS.begin(); it != qTS.end(); ++it)
		std::cout << std::setw(6) << std::setfill('0') << std::setprecision(3) << it->first << "\t" << std::setprecision(5) << it->second << std::endl;
	qTS.clear();

	return true;
}
