#pragma once
#include "Globals.hpp"
#include "ModelGlobals.hpp"


bool InitializePrecipitationModule(ModelParameters const & params);
Vector_f64 ComputePreciptationVector(double time, ModelParameters const & params);