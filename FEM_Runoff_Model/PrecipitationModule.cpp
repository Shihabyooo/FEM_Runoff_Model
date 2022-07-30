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

//the main maps above must be incremented/updated every time step, however, the functions of this class are called several times\
in the iterative solution loop, we store the changes of each call in a temporary map, once the model finishes with iterative solution\
and before it proceeds to the next step, we apply (add or assign) the temp maps to the main maps.
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

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
		if (!InitElementEntries(it->second.id))
			LogMan::Log("WARNING! Failed to create precipitation/moisture tracking entry for element " + std::to_string(it->second.id), LOG_WARN);
	
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
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		//This one is not cummulative.
		currentPrecipitationRate.at(it->second.id) = tempCurrentPrecipRate.at(it->second.id);

		//This one's accumulation is handled by functions that use it.
		totalDirectRunoff.at(it->second.id) = tempTotalDirectRunoff.at(it->second.id);

		totalPrecipitation.at(it->second.id) += tempTotalPrecipitation.at(it->second.id);
		totalLoss.at(it->second.id) += tempTotalLoss.at(it->second.id);
		soilMC.at(it->second.id) += tempSoilMC.at(it->second.id);
	}
}

//Precipitation returned as meters per hour. This is the "total" preciptation , before loss is subtracted.
//current impl doesn't need element, but later it would.
double GetPrecipitationRateAtTime(double time, ModelParameters const & params, Triangle const & element)
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
double GetCummulativePrecipitationAtTime(double time, ModelParameters const & params, Triangle const & element)
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

//return Pe rate in m/hr.
double ComputeInitConstPe(double precipRate, double timeStep, InitialAndConstantParams const * params, Triangle const & element)
{
	double totalPrecip = precipRate * timeStep * 1000.0; //convert back to mm

	//initial loss
	double mcDeficit = Max(params->initialLoss - soilMC.at(element.id), 0.0);
	double remainingPrecip = Max(totalPrecip - mcDeficit, 0.0);

	tempSoilMC.at(element.id) = soilMC.at(element.id) + totalPrecip - remainingPrecip;

	//recompute rate, subtract const rate and return result
	return Max((remainingPrecip / timeStep) - params->constRate, 0.0) / 1000.0;
}

//return Pe rate in m/hr.
double ComputeSCSCNPe(double time, ModelParameters const & params, Triangle const & element)
{
	//Pe = (P - 0.2 S)^2 / ( P + 0.8 S)
	//S = (25400 - 254 CN) / CN
	//This is only carried out if P > initial loss.

	double totalPrecip = GetPrecipitationRateAtTime(time, params, element) * params.timeStep * 1000.0; //TODO total precip rate already computed in calling function. Why compute again?

	double s = (25400.0 - 254.0 * static_cast<double>(params.scsCN)) / static_cast<double>(params.scsCN);
	double initialLoss = 0.2 * s;

	double mcDeficit = Max(initialLoss - soilMC.at(element.id), 0.0);
	double remainingPrecip = Max(totalPrecip - mcDeficit, 0.0);
	tempSoilMC.at(element.id) = soilMC.at(element.id) + totalPrecip - remainingPrecip;

	if (remainingPrecip <= 0.0)
		return 0.0;

	double cummulativePrecip = GetCummulativePrecipitationAtTime(time, params, element) * 1000.0;
	double pe = pow(cummulativePrecip - 0.2 * s, 2.0) / (cummulativePrecip + 0.8 * s);
	tempTotalDirectRunoff.at(element.id) = pe;

	//the incremental precipitation = Pe computed above - Pe up to last timeStep.
	//The second term is stored in the totalDirectRunoff map.
	double rate = Max(pe - totalDirectRunoff.at(element.id), 0.0) / params.timeStep / 1000.0;

	return rate;
}

double ComputeEffectivePreciptationAtTime(double time, ModelParameters const & params, Triangle const & element)
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
	tempTotalPrecipitation.at(element.id) = currentPrecipRate * element.Area() * params.timeStep;
	tempTotalLoss.at(element.id) = (currentPrecipRate - effectiveRate) * element.Area() * params.timeStep;
	
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

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		double newPrecipitation = ComputeEffectivePreciptationAtTime(time + params.timeStep, params, it->second);
		double elementContrib = (it->second.Area() / 3.0)  * ((1.0 - params.femOmega) * currentPrecipitationRate[it->second.id] + params.femOmega * newPrecipitation);
		elementContrib *= params.timeStep;

		tempCurrentPrecipRate.at(it->second.id) = newPrecipitation;

		//Add element's contribution to the global vector.
		for (int intNodeID = 0; intNodeID < 3; intNodeID++)
			result[it->second.VertexID(intNodeID)] += elementContrib;
	}

	return result;
}
