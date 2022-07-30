#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

#include "Globals.hpp"

//TODO CSV parsing should be unified, instead of having several functions with similar code.

//Must be defined in FileIO not in Global because memory freeing requires STB functions
struct Image
{
public:
	Image();
	Image(std::string const & path);
	~Image();

	bool Load(std::string const & path);
	void Unload();

	unsigned char * bitmap = NULL;
	int width = 0, height = 0, samples = 0;
};

namespace FileIO
{
	bool FileExists(std::string const & path);

	bool LoadTimeSeries(std::string const & path, TimeSeries & timeSeries, unsigned int firstLinesToSkip = 0);

	bool LoadCSV(std::string const & path, std::vector<double> & output, unsigned int firstLinesToSkip = 0);
	bool LoadCoordinatePairsCSV(std::string const & path, std::vector<Vector2D> & output, unsigned int firstLinesToSkip = 0);

	bool LoadRaster(std::string const & path, int * outRasterID, Matrix_f64 const ** outBitmapPtr); 
	void UnloadRaster(int rasterID);
	
	bool GetRasterMappingParameters(int rasterID,
									Vector2Int & outDimensions,
									int & outSamples,
									bool & outIsUTM,
									double *** outTiePoints, //Will allocate double[2][3] and point this pointer to it.
									double ** outPixelScale);//will allocate double[3] and point this pointer to it.

	Image LoadImage(std::string const & path);

	bool LoadVectorPath(std::string const & path, std::vector<Vector2D> & output);

	bool InitLogFile();
	void CloseLogFile();
	bool WriteToLog(LogEntry const & newEntry);

	bool InitOutputFile(std::string const & modelName);
	void CloseOutputFile();
	bool WriteOutputFrame(double time, Vector_f64 const & heads, Vector_f64 const & qX, Vector_f64 const & qY);
}