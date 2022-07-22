#pragma once
#include "PrecipitationModule.hpp"
#include "SpatialDataModule.hpp"


std::unordered_map<size_t, double> soilMC; //the current soil moisture content, defined per element (the size_t = element id), and\
											//the double is the MC in mm. For -e.g.- InitConst, the double will start from zero\
											//and grow until initial loss is satisfied.

void ClearData()
{
	soilMC.clear();
}

bool InitializePrecipitationModule(ModelParameters const & params)
{
	LogMan::Log("Initializing precipitation module");

	void ClearData();

	switch (activeMeshType)
	{
	case ElementType::triangle:
		if (params.lossModel != LossModel::none)
			for (auto it = triangles.begin(); it != triangles.end(); ++it)
				if (!soilMC.insert({ it->second.id, 0.0 }).second)
					LogMan::Log("WARNING! Failed to create moisture content tracking entry for element " + std::to_string(it->second.id), LOG_WARN);
		break;
	case ElementType::rectangle:
		if (params.lossModel != LossModel::none)
			for (auto it = rectangles.begin(); it != rectangles	.end(); ++it)
				if (!soilMC.insert({ it->second.id, 0.0 }).second)
					LogMan::Log("WARNING! Failed to create moisture content tracking entry for element " + std::to_string(it->second.id), LOG_WARN);
		break;
	default: //shouldn't happen
		LogMan::Log("ERROR! Undefined mesh type in InitializePrecipitationModule()", LOG_ERROR);
		return false;
	}
	
	return true;
}


double ComputeInitConstPe(double precipRate, double timeStep, InitialAndConstantParams const * params, Element const * element)
{
	double totalPrecip = precipRate * timeStep;

	//initial loss
	double mcDeficit = Max(params->initialLoss - soilMC[element->id], 0.0);
	double remainingPrecip = Max(totalPrecip - mcDeficit, 0.0);
	soilMC[element->id] += totalPrecip - remainingPrecip; 
	
	//recompute rate, subtract const rate and return result
	return Max((remainingPrecip / timeStep) - params->constRate, 0.0);
}


//Precipitation returned as meters per hour. This is the "total" preciptation , before loss is subtracted.
//current impl doesn't need element, but later it would.
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

double ComputeEffectivePreciptation(double time, ModelParameters const & params, Element const * element)
{
	double currentPrecipRate = GetCurrentPrecipitation(time, params, element);

	switch (params.lossModel)
	{
	case LossModel::none:
		return currentPrecipRate;
	case LossModel::initialConst:
		return ComputeInitConstPe(currentPrecipRate, params.timeStep, (InitialAndConstantParams *)params.lossModelParams, element);
	default: //shouldn't happen
		return 0.0;
	}
}

Vector_f64 ComputePreciptationVector(double time, ModelParameters const & params)
{
	//Preciptation Vector is the last term of the RHS of the formulation. i.e. dT * {Beta} * ((1 - omega) * Pe_t + omega * Pe_t+dt)

	//Beta matrix for each element is a 3x1 vector for triangles (4x1 for rectangles)
	//{Beta_e} = A/3 *	|	1	|
	//					|	1	|
	//					|	1	|
	//Where A is the area of element.

	//init outVector to zeroes before doing anything.
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
