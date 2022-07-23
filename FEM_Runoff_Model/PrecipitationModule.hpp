#pragma once
#include "Globals.hpp"
#include "ModelGlobals.hpp"


bool InitializePrecipitationModule(ModelParameters const & params);
double GetWatershedCumulativePrecipitationVolume();
double GetElementCumulativePrecipitationVolume(size_t elementID);

double GetWatershedCumulativeLossVolume();
double GetElementCumulativeLossVolume(size_t elementID);

void AppendLastTimeStepPrecipitationVariables();

Vector_f64 ComputePreciptationVector(double time, ModelParameters const & params);