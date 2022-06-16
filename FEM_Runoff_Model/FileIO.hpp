#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Globals.hpp"

bool LoadCSV(std::string const & path, std::vector<float> & output, unsigned int const & firstLinesToSkip = 0);
bool LoadCoordinatePairsCSV(std::string const & path, std::vector<Vector2> & output, unsigned int const & firstLinesToSkip = 0);

bool LoadRaster(std::string const & path, void * output); //output will be a pointer to a matrixPP_f32