#pragma once
#include "Globals.hpp"
#include "ModelGlobals.hpp"

extern Vector_f64 nodeSlope, nodeManning, nodeFDR;

bool LoadInputRasters(ModelParameters const & params);
void UnloadAllRasters();
double FDR2Angle(int fdr);
bool CacheManningCoefficients(ModelParameters const & params);
bool CacheSlopes(ModelParameters const & params);