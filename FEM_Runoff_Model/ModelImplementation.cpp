#pragma once
#include "ModelImplementation.hpp"
#include "PrecipitationModule.hpp"

std::unordered_map<size_t, Triangle> triangles;
std::vector<Vector2D> nodes;
std::vector<size_t> boundaryNodes;
Vector2D nodesSW, nodesNE;
Vector2D shedSW, shedNE;
std::vector<Vector2D> shedBoundary;

Vector_f64 heads;
std::vector<std::pair<double, double>> outputTimeSeries; //Time series for resulting head at outlet

//TODO have to vector_f64 to store computed nodal velocities, and a function to compute the newheads nodal velocities every \
iteration. newHead velos are then moved to oldHead velos after internal iterations are done.
std::pair<double, double> ComputeVelocityComponents(size_t nodeID, Vector_f64 waterElevation) //waterElevation = nodeElevation + head
{
	double head = waterElevation[nodeID];
	double velocity = 3600.0 * sqrt(nodeSlope[nodeID]) * pow(head, 2.0 / 3.0) / nodeManning[nodeID];
	double angle = FDR2Angle(lround(nodeFDR[nodeID]));

	double u = velocity * cos(angle);
	double v = velocity * sin(angle);

	return std::pair<double, double>(u, v);
}

Matrix_f64 ElementKMatrix(Triangle const & tri, std::pair<double, double> uv[3])
{
	//						|	ui (yj - yk) + vi (xk - xj)		ui (yk - yi) + vi (xi - xk)		ui (yi - yj) + vi (xj - xi)	|
	//[K_element] = 1/6	*	|	uj (yj - yk) + vj (xk - xj)		uj (yk - yi) + vj (xi - xk)		uj (yi - yj) + vj (xj - xi)	|
	//						|	uk (yj - yk) + vk (xk - xj)		uk (yk - yi) + vk (xi - xk)		uk (yi - yj) + vk (xj - xi)	|

	//The 1/6 scalar is mulitplied in the calling function.

	Vector2D const & i = tri.Node(0);
	Vector2D const & j = tri.Node(1);
	Vector2D const & k = tri.Node(2);

	Matrix_f64 kMat(3, 3);

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

Matrix_f64 ElementCMatrixLumped(int size, double areaDivision)
{
	//Lumped Capcitance Matrix for both element type defined as a diagonal matrix with diag value equal to element area / number of nodes\
	i.e. "areaDivision."

	Matrix_f64 cMat(size, size);

	for (int i = 0; i < size; i++)
		cMat[i][i] = areaDivision;

	return cMat;
}

//Returns either consistent or lumped form, depending on setting in params.
Matrix_f64 ElementCMatrix(ModelParameters const & params, Triangle const & tri)
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

Matrix_f64 ComputeGlobalCoefficientsMatrix(ModelParameters const & params, Vector_f64 const & newHeads)
{
	//[C] + w*dt*[K]
	Matrix_f64 coefMat(nodes.size(), nodes.size());
	Matrix_f64 kMat, cMat;

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		std::pair<double, double> uv[3] = { ComputeVelocityComponents(it->second.VertexID(0), newHeads),
											ComputeVelocityComponents(it->second.VertexID(1), newHeads),
											ComputeVelocityComponents(it->second.VertexID(2), newHeads) };

		kMat = ElementKMatrix(it->second, uv) * (params.timeStep * params.femOmega / 6.0);
		cMat = ElementCMatrix(params, it->second);

		for (int row = 0; row < 3; row++)
			for (int column = 0; column < 3; column++)
				coefMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] + kMat[row][column]);
	}

	return coefMat;
}

Matrix_f64 ComputeGlobalConductanceMatrix(ModelParameters const & params)
{
	//[C] - dt*(1-w)*[K]
	Matrix_f64 condMat(nodes.size(), nodes.size());
	Matrix_f64 kMat, cMat;

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		std::pair<double, double> uv[3] = { ComputeVelocityComponents(it->second.VertexID(0), heads),
											ComputeVelocityComponents(it->second.VertexID(1), heads),
											ComputeVelocityComponents(it->second.VertexID(2), heads) };

		kMat = ElementKMatrix(it->second, uv) * (params.timeStep * (1.0 - params.femOmega) / 6.0);
		cMat = ElementCMatrix(params, it->second);

		for (int row = 0; row < 3; row++)
			for (int column = 0; column < 3; column++)
				condMat[it->second.VertexID(row)][it->second.VertexID(column)] += (cMat[row][column] - kMat[row][column]);
	}

	return condMat;
}

//[A]{x} = {b}
void AdjustForBoundaryConditions(Matrix_f64 & aMat, Vector_f64 & bVec)
{
	////Using https://finite-element.github.io/7_boundary_conditions.html
	////Or Method 2 in Rao 2018.
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
	{
		bVec[*it] = 0.0;
		for (size_t i = 0; i < aMat.Columns(); i++)
			aMat[*it][i] = 0.0;
		aMat[*it][*it] = 1.0;
	}
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

	if (!InitializePrecipitationModule(params))
		return false;

	//Special consideration. Since the boundary nodes list from the mesh generator includes our exit node, we have to\
	 manually remove it.
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

	outputTimeSeries.clear();
#pragma endregion

#pragma region Test
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

	std::cout << "\n===================================================\n";
	std::cout << "Elevations, Slopes, Manning roughness coef";
	std::cout << "\n===================================================\n";
	std::cout << "node |   n  |  Slope  |  FDR  \n";
	for (size_t i = 0; i < nodes.size(); i++)
		std::cout << std::fixed << std::setw(4) << std::setprecision(4) << i << " : " <<
					nodeManning[i] << " | " << nodeSlope[i] << " | " << nodeFDR[i]<< std::endl;
#pragma endregion

	//Allocate the heads vector and initialize it to zero (dry bed initial conditions).
	heads = Vector_f64(nodes.size());

	//Loop from start time to end time
	double time = params.startTime;

	//to compute some statistics about error.
	//Finish implementing error statistics.
	std::vector<double> internalResidualsLog;
	std::vector<double> solverResidualsLog;

	LogMan::Log("Starting simulation loop");

	while (time <= params.endTime)
	{
		std::cout << "\n\n===================================================\n";
		LogMan::Log("At T= " + std::to_string(time));
		std::cout << "===================================================\n";

		//Initialize and assume values for newheads. 
		Vector_f64 newHeads(nodes.size());
		//To speed up reaching a solution, only assume non-zero heads at internal nodes.
		for (size_t i = 0; i < nodes.size(); i++)
			if (!IsBoundaryNode(i))
				newHeads[i] = heads[i] + 0.001;

		//external loop
		for (size_t i = 0; i <= params.maxExternalIterations; i++)
		{
			//RHS = [GlobalConductanceMat] * {h_0} + precipComponent
			Vector_f64 RHS = ComputeGlobalConductanceMatrix(params) * heads + ComputePreciptationVector(time, params);;
			Matrix_f64 coeffMat = ComputeGlobalCoefficientsMatrix(params, newHeads);
			AdjustForBoundaryConditions(coeffMat, RHS);

			Vector_f64 residuals; //needed by solver
			Vector_f64 fixedNewH;

			if (!Solve(coeffMat, RHS, fixedNewH, residuals, params))
			{
				LogMan::Log("ERROR! Internal solver error.", LOG_ERROR);
				return false;
			}

			//force computed newHeads to be positive (negative values will break the model, as pow(negative number) = NaN)
			for (size_t j = 0; j < fixedNewH.Rows(); j++)
				fixedNewH[j] = Max(fixedNewH[j], 0.0);

			solverResidualsLog.push_back(residuals.SumAbs());
			internalResidualsLog.push_back((newHeads - fixedNewH).SumAbs());

			if (internalResidualsLog.back() <= params.externalResidualTreshold) //Good result, we can break the loop now.
			{
				newHeads = fixedNewH;
				break;
			}
			else if (i >= params.maxExternalIterations) //While not always a bad thing (e.g. if threshold was unrealistically small), we should warn user regardless.
				LogMan::Log("Reached maxInternalIterations without reaching appropriate h", LOG_WARN);

			newHeads = fixedNewH;
		}

		AppendLastTimeStepPrecipitationVariables(); //See PrecipitationModule.
		
		double q = sqrt(nodeSlope[params.outletNode]) * pow(heads[params.outletNode], 5.0 / 3.0) / nodeManning[params.outletNode];
		std::string timeStepResults = "Head = " + std::to_string(heads[params.outletNode]) + "m - q = " + std::to_string(q) + "m2/s";
		LogMan::Log(timeStepResults);
		outputTimeSeries.push_back(std::pair<double, double>(time, heads[params.outletNode]));

		heads = newHeads;
		time += params.timeStep;
	}

	LogMan::Log("Finished simulation!", LOG_SUCCESS);

	double totalInput = GetWatershedCumulativePrecipitationVolume(); //cubic meters.
	double totalLoss = GetWatershedCumulativeLossVolume(); //cubic meters.
	double totalOutput = 0.0; //cubic meters.
	double flowWidth = 4.0 *  abs((triangles[0].Centroid() - triangles[0].Node(0)).x); //TODO improve flow width computations. Currently only valid for synthetic watershed.

	std::cout << " time \t Q (cms)\n";
	for (auto it = outputTimeSeries.begin(); it != outputTimeSeries.end(); ++it)
	{
		double q = sqrt(nodeSlope[params.outletNode]) * pow(it->second, 5.0 / 3.0) / nodeManning[params.outletNode];
		double Q = q * flowWidth;
		std::cout << std::setw(6) << std::setfill('0') << std::setprecision(3) << it->first << "\t" << std::setprecision(5) << Q << std::endl;
		totalOutput += Q * params.timeStep * 3600.0;
	}
	
	LogMan::Log("Total precipitation: " + std::to_string(totalInput) + " cubic meters.");
	LogMan::Log("Total loss: " + std::to_string(totalLoss) + " cubic meters.");
	LogMan::Log("Total runoff: " + std::to_string(totalOutput) + " cubic meters.");

	return true;
}