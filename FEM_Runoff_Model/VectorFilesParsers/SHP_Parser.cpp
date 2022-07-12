#include "SHP_Parser.h"

SHPParser::SHPParser()
{
	isPathLoaded = false;
	verts = NULL;
	vertsCount = NULL;
}

SHPParser::~SHPParser()
{
	UnLoadGeometry();
}

bool SHPParser::LoadGeometry(std::string const & fileName)
{
	std::string fileNamePrefix = RemoveFileExtension(fileName);

	if (!CheckFileExistance(fileNamePrefix))
	{
		std::cout << "ERROR! One or more of the required ShapeFile files is missing.";
		return false;
	}

	std::cout << "Attempting to load parameters from SHX file." << std::endl; //test
	if (!LoadSHPParameters(fileNamePrefix))
	{
		std::cout << "ERROR! Something went wrong while reading shapefile parameters from .SHX file." << std::endl; //TODO enrich this check
		return false;
	}
	std::cout << "Succesfully loaded parameters from SHX file." << std::endl; //test

	AllocateVertsArray();

	std::cout << "Attempting to extract paths from SHP file." << std::endl; //test
	if (!ExtractPaths(fileNamePrefix))
	{
		std::cout << "ERROR! Something went wrong while extracting the paths" << std::endl; //TODO enrich this check
		return false;
	}
	
	isPathLoaded = true;
	std::cout << "Succesfully loaded paths from SHP file." << std::endl; //test
	
	
	
	std::cout << "Attempting to load CRS info from .prj file." << std::endl; //test	
	if(!LoadCRSDetails(fileNamePrefix)) //Note: Lack if CRS details, while the could be problematic for the running of the ProfileMaker, is  not program breaking and could be supplied/estimated by other means.
	{
		std::cout << "Failed to load CRS info from .prj file." << std::endl; //test
	}
	std::cout << "Successfully loaded CRS info from .prj file." << std::endl; //test

	return true;
}

CRS SHPParser::GeometryCRS()
{
	return geometryCRS;
}

unsigned int SHPParser::UTMZone()
{
	if (geometryCRS != CRS::UTM)
	{
		std::cout << "WARNING! Attempting to access UTM zone details for a non-UTM geoemtry." << std::endl;
		return 0;
	}
	
	return zone;
}

bool SHPParser::IsNorthernHemisphere()
{
	if (geometryCRS != CRS::UTM)
	{
		std::cout << "WARNING! Attempting to access UTM zone details for a non-UTM geoemtry." << std::endl;
		return true;
	}

	return isNorthernHemisphere;
}

void SHPParser::UnLoadGeometry()
{
	if (verts != NULL) //Doesn't delete already check for NULLity?
	{
		for (int i = 0; i < pathsCount; i++)
			verts[i].~Matrix_f64();

		delete[] verts;
		verts = NULL;
	}

	if (vertsCount != NULL)
	{
		delete[] vertsCount;
		vertsCount = NULL;
	}
	
	pathsCount = 0;
	isPathLoaded = false;
}

const std::string SHPParser::RemoveFileExtension(const std::string &fileName) const
{
	return fileName.substr(0, fileName.length() - 4);
}

bool SHPParser::CheckFileExistance(const std::string &fileNamePrefix) const
{
	//In this method, we will check for the existence of all of our needed files. For now, we need the SHX and SHP files.
	//We'll have each failure set a more global bool to false instead of immediatly returning to simulate an error log message (i.e. so the user can know which files are missing).
	std::ifstream testStream;
	bool state = true;
	std::string shpFileComponent = fileNamePrefix + ".shp";
	std::string shxFileComponent = fileNamePrefix + ".shx";
	
	testStream.open(shpFileComponent, std::ios::binary | std::ios::in);
	if (!testStream.is_open())
	{
		std::cout << "ERROR! SHP file is missing." << std::endl;
		state = false;
	}
	testStream.close();
	
	testStream.open(shxFileComponent, std::ios::binary | std::ios::in);
	if (!testStream.is_open())
	{
		std::cout << "ERROR! SHX file is missing." << std::endl;
		state = false;
	}

	testStream.close();


	return state;
}

void SHPParser::AllocateVertsArray()
{
	verts = new Matrix_f64[pathsCount];

	for (int i = 0; i < pathsCount; i++)
		verts[i] = Matrix_f64(vertsCount[i], 2);
	
}

bool SHPParser::LoadSHPParameters(const std::string &fileNamePrefix)
{
	//This class is a single-pass-one, meaning we only scane the file once, load content into memory and be done with it. So we have no use for the SHX file's intended use of "fast seeking."
	//We will only use it figure out how many shapes the SHP contains, and the vertex count of each shape

	std::string shxFileName = fileNamePrefix + ".shx";
	std::ifstream shxFile;
	std::cout << "Attempting to open " << shxFileName.c_str() << "\n";
	shxFile.open(shxFileName, std::ios::binary);

	if (!shxFile.is_open())//redundant, we alread checked that we could open the file in CheckFileExistence(), but still...
	{
		std::cout << "Error: Could not open " << shxFileName.c_str() << "!\n\n";
		shxFile.close(); //necessary?
		return false;
	}

	char byte[4]; //a 4 byte buffer to hold a single int32, will be used several times bellow.

	//We now check that our Shape Type is a polyline by checking that bytes 32-35 of the header are equal to 3 (the code for polylines).
	shxFile.seekg(32, shxFile.beg);
	shxFile.read(byte, sizeof(byte));
	//std::cout << "The type code of the contained geometry is: " << BytesToInt32(byte, false) << std::endl; //test

	if (BytesToInt32(byte, false) != 3)
	{
		std::cout << "ERROR! The provided SHP does not contain a polyline geometry." << std::endl;
		shxFile.close();
		return false;
	}

	//Calculate the number of geometries in the SHP:
	shxFile.seekg(24, shxFile.beg); //go to location of file length data in header
	shxFile.read(byte, sizeof(byte)); //Read the shx file length

	//The shx file length is measured in WORDs, the file length include the header's size (fixed 100 bytes). Each geometry record afterwards is at a fixed 8 bytes (2 x 4bytes int32).
	//i.e. the number of geometries in the shape file = ((2 * file length) - 100) / 8

	pathsCount = ((2 * BytesToInt32(byte, true)) - 100) / 8; //in theory, the result of the outer brackets should always be devisible by 16, giving perfect ints...
	std::cout << "No. of records found: " << pathsCount << std::endl; //test

	if (pathsCount < 1) //if there are no shapes in the SHP, there is no point in continuting this process.
	{
		std::cout << "ERROR! No geometries were found in the provided SHP file." << std::endl;
		shxFile.close();
		return false;
	}

	//now that we know how many geometries our SHP has, we can allocate our vertsCount array.
	vertsCount = new long int[pathsCount];

	//Now to check the metadata of the actual geometries
	shxFile.seekg(100, shxFile.beg); //the header is fixed at 100 bytes, we simply skip it since we don't have a use for anything other than the Shape Type and File Length


	for (int i = 0; i < pathsCount; i++)
	{
		shxFile.read(byte, sizeof(byte)); //the first 4 bytes of a record header contains its offset, no need for them now.
		//TODO replace the read above with a seekg(sizeof(byte), shxFile.cur) (or std::ios::cur if that doesn't work).

		//Assuming each record has only one part (a requirement of this program), the coordinates of the record start 48 bytes from the record data start location.
		//The length of the record is measured in WORDs (2 bytes), the coordiantes are in xy pairs stored as doubles (i.e. each pair is 2 * 8bytes in size).
		//e.g. for a record of length = 104. Length in bytes = 104*2 = 208 bytes. Coords size = 208 - 48 = 160 bytes. Coords pair count = 160 / (8 * 2) = 10 coords.
		//In other words, the no.of vertices = no. of coord pairs = ((2 * record length) - 48) / 16.

		shxFile.read(byte, sizeof(byte)); //the second 4 bytes of a record header contains its length, we use this to calculate how many vertices a shape has.

		long int noOfVerts = ((2 * BytesToInt32(byte, true)) - 48) / 16; //Again: in theory, the result of the outer brackets should always be devisible by 16, giving perfect ints...
		vertsCount[i] = noOfVerts;

		if (noOfVerts < 2)
		{
			std::cout << "Error in verts count for shape no. " << i << ", vertsCount: " << vertsCount[i] << std::endl; //test
			delete[] vertsCount;
			vertsCount = NULL;
			shxFile.close();
			return false;
		}
	}
	
	for (int i = 0; i < pathsCount; i++) //test
		std::cout << "Shape no. " << i << ", vertsCount: " << vertsCount[i] << std::endl; //test

	shxFile.close();
	return true;
}

bool SHPParser::LoadCRSDetails(const std::string & fileNamePrefix)
{
	std::string prjFileName = fileNamePrefix + ".prj";
	std::ifstream prjFile;
	
	std::cout << "Attempting to open " << prjFileName.c_str() << "\n";
	prjFile.open(prjFileName);

	if (!prjFile.is_open())
	{
		std::cout << "ERROR! Could not open " << prjFileName.c_str() << "!\n\n";
		prjFile.close(); //necessary?
		return false;
	}

	char crsType[6]; //According to the WKT standard, the file should begin with a 6 character, capitalized string indicating the CRS type: PROJCS, GEOGCS or GEOCCS for projected, geoegraphd and geocentric CRS respectively.
	prjFile.read(crsType, sizeof(crsType));

	std::string _crsType(crsType); //easier to compare strings than char*

	if (_crsType == "PROJCS")
	{
		std::cout << "Loaded SHP uses a projected CRS" << std::endl; //test

		prjFile.read(crsType, 2); //the next two chars are [ and ", we can safely discard them

		//Extract projection name
		char charBuffer = ' ';
		std::string crsName = "";
		while (charBuffer != '"')
		{
			prjFile.read(&charBuffer, sizeof(charBuffer));

			crsName += charBuffer;

			if (prjFile.eof())
			{
				std::cout << "ERROR! Reached EoF before parsing a correct CRS name.";
				geometryCRS = CRS::undefined;
				return false;
			}
		}

		//check that projection is UTM.
		if (crsName.length() < 12 && crsName.substr(0, 12) != "WGS_1984_UTM" && crsName.substr(0, 12) != "WGS 84 / UTM") //First test to ensure following two don't read beyond range,
		{																													//second is ESRI WKT, third is OGC WKT
			std::cout << "ERROR! Projection used is not WGS84 UTM. Projection used is: " << crsName.substr(0, 12).c_str() << std::endl;
			geometryCRS = CRS::undefined;
			return false; //though technically, this is a not a failure in "parsing" the CRS.
		}

		//Extract Zone
		//Remember that crsName includes the delimiter quotation mark as well.
		//Also note that the WKT names are annoying in that they don't pad the single digit zone number with zero, so the zone segment of the name could be 3 chars or 2 chars.
		char _hemisphere = crsName[crsName.length() - 2];
		if (_hemisphere == 'N')
			isNorthernHemisphere = true;
		else if (_hemisphere == 'S')
			isNorthernHemisphere = false;
		else
		{
			std::cout << "ERROR! Could not determine the hemisphere of the CRS's UTM zone.";
			geometryCRS = CRS::undefined;
			return false;
		}
		
		std::string _zoneString = crsName.substr(crsName.length() - 4, 2);
		if (_zoneString[0] == ' ' || _zoneString[0] == '_')
			_zoneString[0] = '0';

		unsigned int _zone = atoi(_zoneString.c_str());

		if (_zone < 1 || _zone > 60)
		{
			std::cout << "ERROR! Could not determine the CRS's UTM zone.";
			geometryCRS = CRS::undefined;
			return false;
		}
		zone = _zone;

		geometryCRS = CRS::UTM;
		std::cout << "Set UTM CRS to Projected UTM, with zone" << zone << (isNorthernHemisphere? "N" : "S") << std::endl;//test
	}
	else if (_crsType == "GEOGCS") 
	{
		std::cout << "Loaded SHP uses a geographic CRS" << std::endl; //test
		//TODO check that the CRS is WGS84 here
		geometryCRS = CRS::WGS84;
	}
	else if (_crsType == "GEOCCS") //GEOGS and GEOCCS will be treated the same here.
	{
		std::cout << "Loaded SHP uses a geocentric CRS" << std::endl; //test
		//TODO check that the CRS is WGS84 here
		geometryCRS = CRS::WGS84;
	}
	else
	{
		geometryCRS = CRS::undefined;
		std::cout << "ERROR! Could not determine the SHP's CRS." << std::endl; //test
	}

	prjFile.close();
	return true;
}

bool SHPParser::ExtractPaths(const std::string &fileNamePrefix)
{
	std::string shpFileName = fileNamePrefix + ".shp";
	std::ifstream shpFile;
	shpFile.open(shpFileName, std::ios::binary);

	if (!shpFile.is_open())
	{
		std::cout << "Error: Could not open " << shpFileName.c_str() << "!\n\n";
		shpFile.close();
		return false;
	}
	std::cout << "Current Loc: " << shpFile.tellg() << std::endl;

	shpFile.seekg(100, shpFile.beg); //skip the file header

	for (int i = 0; i < pathsCount; i++)
	{
		//std::cout << "Current Loc: " << shpFile.tellg() << std::endl;

		//At the beining of this loop, check that the parts number of the record (in32 starting from 36th Byte) is equal to 1.
		//Skip the record header (which we assume fixed at 48 bytes for a single part geometry, though we could estimate the jump by 44 + 4 x partNo).
		//Now start another loop with vertsCount[i] as limit
			//From here on, the data will be read in double pairs (x and y).
			//First double read will be assigned to verts[i][j][0];
			//Second double read will be assigned to verts[i][j][1];

		shpFile.seekg(8, std::ios::cur); //skip the record's header (fixed 8 bytes, or 2xint32)
		char byte[4]; //honestly, this should be DWORD, or quadBytes, or even plural bytes, not just "byte."
		long int noOfParts;
		
		shpFile.seekg(36, std::ios::cur); //now we are at the NumParts (int32) position of the record's header
		shpFile.read(byte, sizeof(byte));
		noOfParts = BytesToInt32(byte, false);

		if (noOfParts > 1)
			std::cout << "WARNING! No. of parts in geometry " << i << "are greater than one." << std::endl;
		
		//We should be at byte of 40 relative to record conent's begining, the first coord starts at location 44 + 4 * noOfParts, we are at byte 40, so we seek 4 + 4 * noOfParts
		shpFile.seekg((4 + 4 * noOfParts), std::ios::cur);

		for (int j = 0; j < vertsCount[i]; j++)
		{
			double x, y;
			shpFile.read((char*)&x, sizeof(x)); 
			shpFile.read((char*)&y, sizeof(y));

			verts[i][j][0] = x;
			verts[i][j][1] = y;
		}

		//verts[i].DisplayArrayInCLI(); //test
	}

	shpFile.close();
	return true;
}

long int SHPParser::BytesToInt32(const char bytes[4], bool isBigEndian) const
{	
	if (isBigEndian)
		return ((unsigned char)(bytes[0] << 24) | (unsigned char)(bytes[1] << 16) | (unsigned char)(bytes[2] << 8) | (unsigned char)bytes[3]);
	else
		return ((unsigned char)(bytes[3] << 24) | (unsigned char)(bytes[2] << 16) | (unsigned char)(bytes[1] << 8) | (unsigned char)bytes[0]);
}