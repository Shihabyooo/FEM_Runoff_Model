#include <GeoTIFF_Parser.h>

#include "FileIO.hpp"
#include "LogManager.hpp"

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
}

bool FileIO::FileExists(std::string const & path)
{
	if (std::filesystem::exists(path))
		return true;
	else
		return false;
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

bool FileIO::LoadCSV(std::string const & path, std::vector<float> & output, unsigned int firstLinesToSkip)
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

bool FileIO::LoadCoordinatePairsCSV(std::string const & path, std::vector<Vector2> & output, unsigned int firstLinesToSkip)
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
	
	float cachedFloat;
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
			output.push_back(Vector2(cachedFloat, atof(sBuffer.c_str())));
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

bool FileIO::LoadRaster(std::string const & path, int * outRasterID) //TODO handle returning pointer.
{
	LogMan::Log("Attempting to load raster file \"" + path + "\"");
	if (!LoadGeoTIFF(path, outRasterID))
	{
		LogMan::Log("Failed to load raster file \"" + path + "\"", LOG_ERROR);
		return false;
	}

	LogMan::Log("Sucessfully loaded raster file!", LOG_SUCCESS);
}

void FileIO::UnloadRaster(int rasterID)
{
	UnloadGeoTIFF(rasterID);
}

std::ofstream logFile;

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
	std::cout << "Writing to log file\n";//test
	if (!logFile.good() || !logFile.is_open())
	{
		LogMan::Log("ERROR! I/O error.", LOG_ERROR);
		return false;
	}

	logFile << newEntry.content << std::endl;
	logFile.flush(); //flush results to disk
	return true;
}
