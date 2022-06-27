#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

#include "Globals.hpp"

namespace FileIO
{
	bool FileExists(std::string const & path);

	bool LoadTimeSeries(std::string const & path, TimeSeries & timeSeries, unsigned int firstLinesToSkip = 0);

	bool LoadCSV(std::string const & path, std::vector<float> & output, unsigned int firstLinesToSkip = 0);
	bool LoadCoordinatePairsCSV(std::string const & path, std::vector<Vector2> & output, unsigned int firstLinesToSkip = 0);

	bool LoadRaster(std::string const & path, int * outRasterID); //output will be a pointer to a matrixPP_f32
	void UnloadRaster(int rasterID);

	bool InitLogFile();
	void CloseLogFile();
	bool WriteToLog(LogEntry const & newEntry);
}