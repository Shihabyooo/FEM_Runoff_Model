#include "KML_Parser.h"

KMLParser::KMLParser()
{
	verts = NULL;
	pathsCount = 0;
	isPathLoaded = false;
}

KMLParser::~KMLParser()
{
	UnLoadGeometry();
}

bool KMLParser::LoadGeometry(std::string const & fileName)
{
	if (!OpenKMLFile(fileName))
	{
		std::cout << "Failed to load KML: " << fileName.c_str() << std::endl;
		return false;
	}

	if (!ExtractPaths())
	{
		std::cout << "Failed to extract paths from KML: " << fileName.c_str() << std::endl;
		CloseKMLFile();
		return false;
	}

	CloseKMLFile();
	isPathLoaded = true;
	return true;
}

void KMLParser::UnLoadGeometry()
{
	if (verts != NULL)
	{
		for (int i = 0; i < pathsCount; i++)
			verts[i].~Matrix_f64();

		delete[] verts;
		verts = NULL;
	}

	pathsCount = 0;	
	kmlFile.clear();
	isPathLoaded = false;
}

CRS KMLParser::GeometryCRS()
{
	return geometryCRS;
}

bool KMLParser::OpenKMLFile(std::string const & fileName)
{
	kmlFile.open(fileName, std::ios::in | std::ios::binary);
	if (!kmlFile.is_open())
	{
		std::cout << "Error: Could not load KML file!\n\n";
		return false;
	}
	return true;
}

void KMLParser::CloseKMLFile()
{
	if(kmlFile.is_open())
		kmlFile.close();
}

bool KMLParser::ExtractPaths()
{
	kmlFile.clear();
	kmlFile.seekg(0, kmlFile.beg);

	std::unique_ptr<char> tag = AdvanceToNextTag();

	//Find how many lines the kml has -> Allocate matrices accordingely.
	while (!kmlFile.eof() && tag != NULL)
	{
		if (CompareTags(tag.get(), KML_POLYLINE_TAG) || CompareTags(tag.get(), KML_POLYGON_TAG))
			pathsCount++;

		tag = AdvanceToNextTag();
	}
	
	if (pathsCount < 1)
	{
		std::cout << "Not paths were found in the kml file. Wrong data type? (Type must be lines)" << std::endl; //test
		//No memory is initialized. Can safely return.
		return false;
	}

	std::cout << "Found " << pathsCount << " paths in KML file." << std::endl; //test
	
	verts = new Matrix_f64[pathsCount];
	pathsNames = new std::string[pathsCount];
	
	kmlFile.clear();
	kmlFile.seekg(0, kmlFile.beg);
	
	tag = AdvanceToNextTag();

	int counter = 0;
	while (!kmlFile.eof() && tag != NULL)
	{
		tag = AdvanceToNextTag();

		if (CompareTags(tag.get(), "<Placemark>"))
		{
			while (!CompareTags(tag.get(), "</Placemark>") && tag != NULL)
			{
				tag = AdvanceToNextTag(); 
				
				if (CompareTags(tag.get(), KML_NAME_TAG))
				{
					KMLElement element = GetCurrentElementValue();
					
					ExtractNameFromKMLElement(&element, counter); //ExtractNameFromKMLElement handles issues by itself. Naming is not critical, no need to break entire program if it had issues.
				}
				else if (CompareTags(tag.get(), KML_COORDINATES_TAG))
				{
					KMLElement element = GetCurrentElementValue();

					if (!ExtractCoordinatesFromKMLElement(&element, counter))
						return false;
				}
			}
			counter++;
		}
	}

	return true;
}

std::unique_ptr<char> KMLParser::AdvanceToNextTag()
{
	if (!kmlFile.is_open())
		return NULL;

	char letter = ' ';
	unsigned long long int tagBeginLoc, tagEndLoc;

	while (letter != '<') //seek begining of tag;
	{
		kmlFile.read(&letter, sizeof(letter));
		if (kmlFile.eof())
			return NULL;
	}

	tagBeginLoc = kmlFile.tellg();
	tagBeginLoc -= 1; //to accomodate the character we've already read.

	letter = ' ';
	while (letter != '>')
	{
		kmlFile.read(&letter, sizeof(letter));
		if (kmlFile.eof())
			return NULL;
	}

	tagEndLoc = kmlFile.tellg();

	std::unique_ptr<char> tag(new char[tagEndLoc - tagBeginLoc]);

	kmlFile.seekg(tagBeginLoc, kmlFile.beg);

	for (int i = 0; i < tagEndLoc - tagBeginLoc; i++)
	{
		kmlFile.read(&letter, sizeof(letter));
		tag.get()[i] = letter;
	}

	return tag;
}

KMLElement KMLParser::GetCurrentElementValue()
{
	if (!kmlFile.is_open())
		return KMLElement();

	char letter = ' ';
	unsigned long long int valBeginLoc, valEndLoc;

	valBeginLoc = kmlFile.tellg();

	letter = ' ';
	while (letter != '<')
	{
		kmlFile.read(&letter, sizeof(letter));
		if (kmlFile.eof())
			return KMLElement();
	}

	valEndLoc = kmlFile.tellg();
	valEndLoc -= 1;

	KMLElement value;

	value.ptrToData = std::unique_ptr<char>(new char[valEndLoc - valBeginLoc]);
	value.dataLength = valEndLoc - valBeginLoc;

	kmlFile.seekg(valBeginLoc, kmlFile.beg);

	for (int i = 0; i < valEndLoc - valBeginLoc; i++)
	{
		kmlFile.read(&letter, sizeof(letter));
		value.ptrToData.get()[i] = letter;
	}

	return value;
}

bool KMLParser::CompareTags(char const * tag1, char const * tag2) const
{
	if (tag1 == NULL)
	{
		//std::cout << "WARNING! Trying to compare a tag of NULL value" << std::endl;
		return false;
	}

	while (*tag1 == *tag2)
	{
		if (*tag1 == '>' && *tag2 == '>')
			return true;

		tag1++;
		tag2++;
	}
	
	return false;
}

void KMLParser::ExtractNameFromKMLElement(KMLElement * element, int pathID)
{
	if (element->dataLength > 0)
	{
		std::string _name(element->ptrToData.get());
		_name[element->dataLength] = '\0'; //strings must be null-terminated
		pathsNames[pathID] = _name;
	}
	else
	{
		pathsNames[pathID] = "UntitledPath_" + pathID;
	}
}

bool KMLParser::ExtractCoordinatesFromKMLElement(KMLElement * element, int pathID)
{
	int commaDelimitersPerCoordinate = 0;

	//first, determine whether coordinates are saved in pairs or trios.
	for (int i = 0; i < element->dataLength; i++)
	{
		if (element->ptrToData.get()[i] == ' ')
		{
			for (int j = i - 1; j >= 0; j--)
			{
				if (element->ptrToData.get()[j] == ',')
					commaDelimitersPerCoordinate++;
			}
			break;
		}
		else if (i == element->dataLength - 1) //reached end of element without finding a space delimiter
			return false;
	}

	if (commaDelimitersPerCoordinate < 1 || commaDelimitersPerCoordinate > 2)
		return false;

	//Now we read our elements to a temporary Vector.
	Matrix_f64 coord = Matrix_f64(3, 1);
	std::vector<Matrix_f64> tempContainer;
	std::string sBuffer;

	for (int i = 0; i < element->dataLength; i++) //loop through entire element
	{
		if (element->ptrToData.get()[i] == '0' && i < element->dataLength - 1)
			i++;

		for (int j = 0; j < commaDelimitersPerCoordinate + 1; j++) //loop through a single pair/trio.
		{
			sBuffer = "";
			while (element->ptrToData.get()[i] != ',' && element->ptrToData.get()[i] != ' ' && i < element->dataLength) 
			{
				sBuffer += element->ptrToData.get()[i];
				i++;
			}
			coord[j][0] = std::stod(sBuffer.c_str());
			i++;
		}
		
		i--; //the last i++ of the j-loop would jump over the space delimiter, resulting in a removal of the largest significant digit of the following x-coordinates.
		tempContainer.push_back(coord);
	}


	//We then move the data from our vector to our Matrix_f64 *;
	verts[pathID] = Matrix_f64(tempContainer.size(), 2);
	int counter = 0;
	for (std::vector<Matrix_f64>::iterator it = tempContainer.begin(); it < tempContainer.end(); ++it)
	{
		verts[pathID][counter][0] = it->GetValue(0, 0);
		verts[pathID][counter][1] = it->GetValue(1, 0);
		counter++;
	}
	
	return true;
}