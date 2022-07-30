#pragma once
#include <GeoTIFF_Parser.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include "FileIO.hpp"
#include "LogManager.hpp"
#include "KML_Parser.h"
#include "SHP_Parser.h"

bool SkipNLines(std::ifstream & file, unsigned int const & linesToSkip)
{
	for (size_t i = 0; i < linesToSkip; i++)
	{
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (file.eof())
		{
			LogMan::Log("Error parsing file. Unexpected endl of file.", LOG_ERROR);
			file.close();
			return false;
		}
		else if (file.fail() || file.bad())
		{
			LogMan::Log("Error parsing file. Internal error or bad file.", LOG_ERROR);
			file.close();
			return false;
		}
	}
	return true;
}

bool FileIO::FileExists(std::string const & path)
{
	return std::filesystem::exists(path);
}

bool FileIO::LoadTimeSeries(std::string const & path, TimeSeries & timeSeries, unsigned int firstLinesToSkip)
{
	LogMan::Log("Attempting to load Time-series csv file \"" + path + "\"");
	std::ifstream file;
	file.open(path);

	if (!file.is_open())
	{
		LogMan::Log("Failed to open file \"" + path + "\"", LOG_ERROR);
		return false;
	}
	
	char cBuffer = ' ';
	std::string sBuffer = "";
	std::vector<std::pair<size_t, double>> tsBuffer;

	//skip firstLinesToSkip
	SkipNLines(file, firstLinesToSkip);

	long cachedLong;
	bool second = false;

	while (!file.eof())
	{
		file.read(&cBuffer, sizeof(cBuffer));

		if (cBuffer == ',')
		{
			cachedLong = atol(sBuffer.c_str());
			sBuffer = "";
			second = true;
		}
		else if ((cBuffer == '\n' || file.eof()) && second)
		{
			tsBuffer.push_back(std::pair<size_t, double>(cachedLong, atof(sBuffer.c_str())));
			sBuffer = "";
			second = false;
		}
		else
			sBuffer += cBuffer;
	}

	file.close();
	timeSeries = TimeSeries(tsBuffer);
	LogMan::Log("Sucessfully loaded Time-series file!", LOG_SUCCESS);
	return true;
}

bool FileIO::LoadCSV(std::string const & path, std::vector<double> & output, unsigned int firstLinesToSkip)
{
	LogMan::Log("Attempting to load csv file \"" + path + "\"");
	std::ifstream file;
	file.open(path);
	
	if (!file.is_open())
	{
		LogMan::Log("Failed to open file \"" + path + "\"", LOG_ERROR);
		return false;
	}

	char cBuffer = ' ';
	std::string sBuffer = "";

	//skip firstLinesToSkip
	SkipNLines(file, firstLinesToSkip);

	while (!file.eof())
	{
		file.read(&cBuffer, sizeof(cBuffer));

		if (cBuffer == ',' || cBuffer == '\n' || file.eof())
		{
			output.push_back(atof(sBuffer.c_str()));
			sBuffer = "";
		}
		else
			sBuffer += cBuffer;
	}

	file.close();
	LogMan::Log("Sucessfully loaded csv file!", LOG_SUCCESS);
	return true;
}

bool FileIO::LoadCoordinatePairsCSV(std::string const & path, std::vector<Vector2D> & output, unsigned int firstLinesToSkip)
{
	LogMan::Log("Attempting to load csv file \"" + path + "\"");
	std::ifstream file;
	file.open(path);

	if (!file.is_open())
	{
		LogMan::Log("Failed to open file \"" + path + "\"", LOG_ERROR);
		return false;
	}

	char cBuffer = ' ';
	std::string sBuffer = "";

	//skip firstLinesToSkip
	SkipNLines(file, firstLinesToSkip);
	
	double cachedFloat;
	bool second = false;

	while (!file.eof())
	{
		file.read(&cBuffer, sizeof(cBuffer));
		
		if (cBuffer == ',')
		{
			cachedFloat = atof(sBuffer.c_str());
			sBuffer = "";
			second = true;
		}
		else if ((cBuffer == '\n' || file.eof()) && second)
		{
			output.push_back(Vector2D(cachedFloat, atof(sBuffer.c_str())));
			sBuffer = "";
			second = false;
		}
		else
			sBuffer += cBuffer;
	}

	file.close();
	LogMan::Log("Sucessfully loaded csv file!", LOG_SUCCESS);
	return true;
}

bool FileIO::LoadRaster(std::string const & path, int * outRasterID, Matrix_f64 const ** outBitmapPtr)
{
	LogMan::Log("Attempting to load raster file \"" + path + "\"");
	if (!LoadGeoTIFF(path, outRasterID))
	{
		LogMan::Log("Failed to load raster file \"" + path + "\"", LOG_ERROR);
		return false;
	}

	*outBitmapPtr = GetBand(*outRasterID, 0);
	LogMan::Log("Sucessfully loaded raster file!", LOG_SUCCESS);
	return true;
}

void FileIO::UnloadRaster(int rasterID)
{
	UnloadGeoTIFF(rasterID);
}

bool FileIO::GetRasterMappingParameters(int rasterID, Vector2Int & outDimensions, int & outSamples, bool & outIsUTM, double *** outTiePoints, double ** outPixelScale)
{
	//There is an assumption here that the A: raster is strictly tie-and-point and, B: it's either Geographic CRS or UTM.
	//Probably should check and enforce rasters to be of those features at loading...

	GeoTIFFDetails const * geoDetails = GetPointerToGeoTIFFDetails(rasterID);
	TIFFDetails const * tiffDetails = GetPointerToTIFFDetails(rasterID);

	if (geoDetails == NULL || tiffDetails == NULL)
	{
		LogMan::Log("ERROR! Failed to get raster tiff/geotiff parameters!", LOG_ERROR);
		return false;
	}

	outDimensions.x = tiffDetails->width;
	outDimensions.y = tiffDetails->height;
	outSamples = tiffDetails->samplesPerPixel;

	outIsUTM = geoDetails->modelType == 1;

	*outTiePoints = new double *[2];
	(*outTiePoints)[0] = new double[3]();
	(*outTiePoints)[1] = new double[3]();

	*outPixelScale = new double[3]();

	for (int i = 0; i < 3; i++)
	{
		(*outTiePoints)[0][i] = geoDetails->tiePoints[0][i];
		(*outTiePoints)[1][i] = geoDetails->tiePoints[1][i];
		(*outPixelScale)[i] = geoDetails->pixelScale[i];
	}

	return true;
}

Image FileIO::LoadImage(std::string const & path)
{
	Image newImage (path);
	return newImage;
}

std::ofstream logFile;

bool FileIO::LoadVectorPath(std::string const & path, std::vector<Vector2D> & output)
{
	//TODO refactor this

	//This is expensive, memory wise. But since, in the current impl, I don't want to maintain neither a FileParser instance and handles (like\
	I do the GeoTiffParser), I simply instantiate a FileParser only in this scope, then copy the data loaded internally to output,\
	then have the internal copy freed with the object's destruction as it goes out of scope.\
	Too much redundancy, but saves LoC in management and spares me the pain of modifiying the parsers...

	std::string extension = path.substr(path.length() - 4, 4);
	FileParser * parser = NULL;

	if (extension == ".kml")
	{
		LogMan::Log("Attempting to read KML file \"" + path + "\"");
		parser = new KMLParser();
		if (!parser->LoadGeometry(path))
		{
			LogMan::Log("ERROR! Could not load from file \"" + path + "\"", LOG_ERROR);
			delete parser;
			return false;
		}
	}
	else if (extension == ".shp")
	{
		LogMan::Log("Attempting to read SHP file \"" + path + "\"");
		parser = new SHPParser();
		if (!parser->LoadGeometry(path))
		{
			LogMan::Log("ERROR! Could not load from file \"" + path + "\"", LOG_ERROR);
			delete parser;
			return false;
		}
	}
	else
	{
		LogMan::Log("ERROR! Unsupported file format.", LOG_ERROR);
		return false;
	}

	//if we reached this point, parser is definitely instantiated and has loaded something
	output.clear();
	Matrix_f64 const * vertices = parser->GetPathByID(0);
	for (size_t i = 0; i < vertices->Rows(); i++)
	{
		Vector2D coords(vertices->GetValue(i, 0), vertices->GetValue(i, 1));
		
		//check if the CRS isn't UTM, convert if so
		if (parser->GeometryCRS() != CRS::UTM)
			coords = ProjectPoint(coords);
		
		output.push_back(coords);
	}

	delete parser;
	LogMan::Log("Successfully loaded file.", LOG_SUCCESS);
	return true;
}

bool FileIO::InitLogFile()
{
	logFile.open(LOG_FILE_PATH, std::ios_base::app); //open LOG_FILE_PATH and always seek to the end before every write
	
	if (!logFile.is_open())
	{
		LogMan::Log("ERROR! Could not open log file " + std::string(LOG_FILE_PATH), LOG_ERROR);
		return false;
	}
		
	return true;
}

void FileIO::CloseLogFile()
{
	if (logFile.is_open())
		logFile.close();
}

bool FileIO::WriteToLog(LogEntry const & newEntry)
{
	if (!logFile.good() || !logFile.is_open())
	{
		LogMan::Log("ERROR! I/O error.", LOG_ERROR);
		return false;
	}

	logFile << newEntry.content << std::endl;
	logFile.flush(); //flush results to disk
	return true;
}

std::ofstream outputFile;

bool FileIO::InitOutputFile(std::string const & modelName)
{
	std::string outputPath = modelName + "_output.csv";
	outputFile.open(outputPath.c_str(), std::ios_base::trunc);
	if (!outputFile.is_open())
	{
		LogMan::Log("ERROR!C could not open or create output file " + outputPath, LOG_ERROR);
		return false;
	}
	return true;
}

void FileIO::CloseOutputFile()
{
	if (outputFile.is_open())
		outputFile.close();
}

bool FileIO::WriteOutputFrame(double time, Vector_f64 const & heads, Vector_f64 const & qX, Vector_f64 const & qY)
{
	if (!outputFile.good() || !outputFile.is_open())
	{
		LogMan::Log("ERROR! I/O error.", LOG_ERROR);
		return false;
	}

	size_t rows = heads.Rows();

	outputFile << ",,,\n";
	outputFile << "Time:," << time << ",,\n";
	outputFile << "NodeID,Head(m),q-x(m2/hr),q-y(m2/hr)\n";
	for (size_t i = 0; i < rows; i++)
		outputFile << i << "," <<
		heads.GetValue(i) << "," <<
		qX.GetValue(i) << "," <<
		qY.GetValue(i) << "\n";

	outputFile.flush();
	return true;
}

#pragma region Image struct defs

Image::Image()
{
}

Image::Image(std::string const & path)
{
	Load(path);
}

Image::~Image()
{
	Unload();
}

bool Image::Load(std::string const & path)
{
	Unload();
	bitmap = stbi_load(path.c_str(), &width, &height, &samples, 0);

	if (bitmap == NULL)
	{
		std::string errorMsg = "Warning! Could not load image \"" + path + "\". -- " + stbi_failure_reason();
		LogMan::Log(errorMsg.c_str(), LOG_WARN);
		return false;
	}

	return true;
}

void Image::Unload()
{
	if (bitmap != NULL)
		stbi_image_free(bitmap);

	bitmap = NULL;
	
	width = 0;
	height = 0;
	samples = 0;
}

#pragma endregion