#pragma once
#include "PrecipitationModule.hpp"
#include "SpatialDataModule.hpp"

std::unordered_map<size_t, double> soilMC; //the current soil moisture content, defined per element (the size_t = element id), and\
											the double is the MC in mm. For -e.g.- InitConst, the double will start from zero\
											and grow until initial loss is satisfied.

std::unordered_map<size_t, double> totalPrecipitation; //The total precipitation input per element over the course of the simulation,\
														in meters head (to be multiplied with the element's area)

std::unordered_map<size_t, double> totalLoss; //The total loss in the system from infiltration, interception, depression, evaporation\
												and/or EVT.

std::unordered_map<size_t, double> totalDirectRunoff; //cummulative depth of effective precipitation

std::unordered_map<size_t, double> currentPrecipitationRate; //the current precipitaiton rate, defined per element (the size_t = element id), and\
															the double is the rate in m / hr. 


//TODO explain 
std::unordered_map<size_t, double> tempSoilMC;
std::unordered_map<size_t, double> tempTotalPrecipitation;
std::unordered_map<size_t, double> tempTotalLoss;
std::unordered_map<size_t, double> tempTotalDirectRunoff;
std::unordered_map<size_t, double> tempCurrentPrecipRate;

void ClearData()
{
	soilMC.clear();
	totalPrecipitation.clear();
	totalLoss.clear();
	totalDirectRunoff.clear();
	currentPrecipitationRate.clear();

	tempSoilMC.clear();
	tempTotalPrecipitation.clear();
	tempTotalLoss.clear();
	tempTotalDirectRunoff.clear();
	tempCurrentPrecipRate.clear();
}

bool InitElementEntries(size_t id, double initialValue = 0.0)
{
	bool status = true;

	status = status && soilMC.insert({ id, 0.0 }).second;
	status = status && totalPrecipitation.insert({ id, 0.0 }).second;
	status = status && currentPrecipitationRate.insert({ id, 0.0 }).second;
	status = status && totalLoss.insert({ id, 0.0 }).second;
	status = status && totalDirectRunoff.insert({ id, 0.0 }).second;

	status = status && tempSoilMC.insert({ id, 0.0 }).second;
	status = status && tempTotalPrecipitation.insert({ id, 0.0 }).second;
	status = status && tempTotalLoss.insert({ id, 0.0 }).second;
	status = status && tempTotalDirectRunoff.insert({ id, 0.0 }).second;
	status = status && tempCurrentPrecipRate.insert({ id, 0.0 }).second;

	return status;
}

bool InitializePrecipitationModule(ModelParameters const & params)
{
	LogMan::Log("Initializing precipitation module");

	ClearData();

	switch (activeMeshType)
	{
	case ElementType::triangle:
		for (auto it = triangles.begin(); it != triangles.end(); ++it)
			if (!InitElementEntries(it->second.id))
				LogMan::Log("WARNING! Failed to create precipitation/moisture tracking entry for element " + std::to_string(it->second.id), LOG_WARN);
		break;
	case ElementType::rectangle:
		for (auto it = rectangles.begin(); it != rectangles	.end(); ++it)
			if (!InitElementEntries(it->second.id))
				LogMan::Log("WARNING! Failed to create precipitation/moisture  content tracking entry for element " + std::to_string(it->second.id), LOG_WARN);
		break;
	default: //shouldn't happen
		LogMan::Log("ERROR! Undefined mesh type in InitializePrecipitationModule()", LOG_ERROR);
		return false;
	}
	
	return true;
}

double SumCumulativeVolume(std::unordered_map<size_t, double> const & targetDepthRecords)
{
	double sum = 0.0;
	for (auto it = targetDepthRecords.begin(); it != targetDepthRecords.end(); ++it)
		sum += it->second;

	return sum;
}

double GetWatershedCumulativePrecipitationVolume()
{
	return SumCumulativeVolume(totalPrecipitation);
}

double GetElementCumulativePrecipitationVolume(size_t elementID)
{
	if (elementID >= totalPrecipitation.size())
	{
		LogMan::Log("WARNING! Attempting to read precipitation volume rectod with an invalid element ID: " + std::to_string(elementID), LOG_WARN);
		return DBL_MIN;
	}
	
	return totalPrecipitation.at(elementID);
}

double GetWatershedCumulativeLossVolume()
{
	return SumCumulativeVolume(totalLoss);
}

double GetElementCumulativeLossVolume(size_t elementID)
{
	if (elementID >= totalPrecipitation.size())
	{
		LogMan::Log("WARNING! Attempting to read loss volume rectod with an invalid element ID: " + std::to_string(elementID), LOG_WARN);
		return DBL_MIN;
	}

	return totalLoss.at(elementID);
}

void AppendLastTimeStepPrecipitationVariables()
{

	size_t elementCount = 0;

	//Can't loop over elements using iterators, so get the number of elements based on activeMeshType, the loop using standard for loops.
	switch (activeMeshType)
	{
	case ElementType::rectangle:
		elementCount = rectangles.size();
		break;
	case ElementType::triangle:
		elementCount = triangles.size();
		break;
	}

	//Note: This assumes the meshing algorithms produced output that is id'ed sequentially from zero to size()-1.
	for (size_t i = 0; i < elementCount; i++)
	{
		//This one is not cummulative.
		currentPrecipitationRate.at(i) = tempCurrentPrecipRate.at(i);
		
		//This one's accumulation is handled by functions that use it.
		totalDirectRunoff.at(i) = tempTotalDirectRunoff.at(i);

		totalPrecipitation.at(i) += tempTotalPrecipitation.at(i);
		totalLoss.at(i) += tempTotalLoss.at(i);
		soilMC.at(i) += tempSoilMC.at(i);
	}
}

//return Pe rate in m/hr.
double ComputeInitConstPe(double precipRate, double timeStep, InitialAndConstantParams const * params, Element const * element)
{
	double totalPrecip = precipRate * timeStep;

	//initial loss
	double mcDeficit = Max(params->initialLoss - soilMC.at(element->id), 0.0);
	double remainingPrecip = Max(totalPrecip - mcDeficit, 0.0);

	//tempSoilMC.at(element->id) =  totalPrecip - remainingPrecip;
	tempSoilMC.at(element->id) = soilMC.at(element->id) + totalPrecip - remainingPrecip;

	//recompute rate, subtract const rate and return result
	return Max((remainingPrecip / timeStep) - params->constRate, 0.0);
}

//return Pe rate in m/hr.
double ComputeSCSCNPe(double time, ModelParameters const & params, Element const * element)
{
	//Pe = (P - 0.2 S)^2 / ( P + 0.8 S)
	//S = (25400 - 254 CN) / CN
	//This is only carried out if P > initial loss

	double totalPrecip = GetPrecipitationRateAtTime(time, params, element) * params.timeStep; //TODO total precip rate already computed in calling function. Why compute again?

	double s = (25400.0 - 254.0 * static_cast<double>(params.scsCN)) / static_cast<double>(params.scsCN);
	double initialLoss = 0.2 * s;

	double mcDeficit = Max(initialLoss - soilMC.at(element->id), 0.0);
	double remainingPrecip = Max(totalPrecip - mcDeficit, 0.0);

	if (remainingPrecip <= 0.0)
		return 0.0;

	double cummulativePrecip = GetCummulativePrecipitationAtTime(time, params, element);
	double pe = (cummulativePrecip - 0.2 * s) / (cummulativePrecip + 0.8 * s);
	tempTotalDirectRunoff.at(element->id) = pe;

	//the incremental precipitation = Pe computed above - Pe up to last timeStep.
	//The second term is stored in the totalDirectRunoff map.
	double rate = Max(pe - totalDirectRunoff.at(element->id), 0.0) / params.timeStep;
	return rate;
}

//Precipitation returned as meters per hour. This is the "total" preciptation , before loss is subtracted.
//current impl doesn't need element, but later it would.
double GetPrecipitationRateAtTime(double time, ModelParameters const & params, Element const * element)
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

//returns cummulative precip in meters
double GetCummulativePrecipitationAtTime(double time, ModelParameters const & params, Element const * element)
{
	if (params.variablePrecipitation)
	{
		LogMan::Log("Warning! Variable precipitation is not yet implemented!", LOG_WARN);
		return 0.0;
	}
	else
	{
		return params.unitTimeSeries.SampleCummulativePreciptation(time) / 1000.0;
	}
}

double ComputeEffectivePreciptationAtTime(double time, ModelParameters const & params, Element const * element)
{
	double currentPrecipRate = GetPrecipitationRateAtTime(time, params, element);
	double effectiveRate = 0.0;

	switch (params.lossModel)
	{
	case LossModel::none:
		effectiveRate = currentPrecipRate;
		break;
	case LossModel::initialConst:
		effectiveRate = ComputeInitConstPe(currentPrecipRate, params.timeStep, (InitialAndConstantParams *)params.lossModelParams, element);
		break;
	case LossModel::scsCN:
		effectiveRate = ComputeSCSCNPe(time, params,  element);
		break;
	default: //shouldn't happen
		break;
	}

	//Update total precip and loss (after converting it back from rate to depth)
	tempTotalPrecipitation.at(element->id) = currentPrecipRate * element->Area() * params.timeStep;
	tempTotalLoss.at(element->id) = (currentPrecipRate - effectiveRate) * element->Area() * params.timeStep;
	
	//test
	/*if (element->id == 0)
		std::cout << "totalPrecip for elem: " << element->id << " = " << totalPrecipitation.at(element->id) << std::endl;*/

	return effectiveRate;
}

Vector_f64 ComputePreciptationVector(double time, ModelParameters const & params)
{
	//Preciptation Vector is the last term of the RHS of the formulation. i.e. dT * {Beta} * ((1 - omega) * Pe_t + omega * Pe_t+dt)

	//Beta vector  for each element is a 3x1 vector for triangles (4x1 for rectangles, times A/4)
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.

	//init outVector to zeroes before doing anything.
	Vector_f64 result(nodes.size());

	size_t elementCount = 0;

	//Can't loop over elements using iterators, so get the number of elements based on activeMeshType, the loop using standard for loops.
	switch (activeMeshType)
	{
	case ElementType::rectangle:
		elementCount = rectangles.size();
		break;
	case ElementType::triangle:
		elementCount = triangles.size();
		break;
	}

	//Note: This assumes the meshing algorithms produced output that is id'ed sequentially from zero to size()-1.
	for (size_t i = 0; i < elementCount; i++)
	{
		Element const * element = NULL;
		double mult = 0.0;
		
		switch (activeMeshType)
		{
		case ElementType::rectangle:
			element = new Element(rectangles[i]); 
			mult = element->Area() / 4.0;
			break;
		case ElementType::triangle:
			element = new Element(triangles[i]);
			mult = element->Area() / 3.0;

			//test
			/*double denom = 0.0;
			for (int i = 0; i < 3; i++)
				if (!IsBoundaryNode(element->VertexID(i)))
					denom += 1.0;
			mult = element->Area() / denom;
			std::cout << "elem: " << element->id << " - area deno: " << denom << std::endl;*/
			break;
		}

		double newPrecipitation = ComputeEffectivePreciptationAtTime(time + params.timeStep, params, element);
		double elementContrib = mult * ((1.0 - params.femOmega) * currentPrecipitationRate[element->id] + params.femOmega * newPrecipitation);
		elementContrib *= params.timeStep;
		
		tempCurrentPrecipRate.at(element->id) = newPrecipitation;
		//tempTotalDirectRunoff.at(element->id) = newPrecipitation * element->Area() * params.timeStep;

		//Add element's contribution to the global vector.
		for (int intNodeID = 0; intNodeID < element->NodeCount(); intNodeID++)
			result[element->VertexID(intNodeID)] += elementContrib;

		delete element;
	}

	return result;
}
