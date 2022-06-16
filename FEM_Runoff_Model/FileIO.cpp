#include <GeoTIFF_Parser.h>

#include "FileIO.hpp"
#include "LogManager.hpp"


bool SkipNLines(std::ifstream & file, unsigned int const & linesToSkip)
{
	for (int i = 0; i < linesToSkip; i++)
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

bool LoadCSV(std::string const & path, std::vector<float> & output, unsigned int const & firstLinesToSkip)
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

bool LoadCoordinatePairsCSV(std::string const & path, std::vector<Vector2> & output, unsigned int const & firstLinesToSkip)
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

bool LoadRaster(std::string const & path, void * output) //TODO handle returning pointer.
{
	LogMan::Log("Attempting to load raster file \"" + path + "\"");
	if (!LoadGeoTIFF(path))
	{
		LogMan::Log("Failed to load raster file \"" + path + "\"", LOG_ERROR);
		return false;
	}

	LogMan::Log("Sucessfully loaded raster file!", LOG_SUCCESS);
}