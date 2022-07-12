#pragma once
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "Vector_File_Parser.h"

#define KML_POLYLINE_TAG "<LineString>"
#define KML_POLYGON_TAG "<LinearRing>" //It's actually something like <Polygon><outerBoundaryIs><LinearRing>, but since I don't want to modify the code too much...
#define KML_NAME_TAG "<name>"
#define KML_COORDINATES_TAG "<coordinates>"

struct KMLElement
{
public:
	KMLElement()
	{
		ptrToData = NULL;
		dataLength = 0;
	}

	std::unique_ptr<char> ptrToData;
	long long int dataLength;
};

class KMLParser : public FileParser
{
public:
	KMLParser();
	~KMLParser();

	bool LoadGeometry(std::string const & fileName);
	void UnLoadGeometry();
	
	virtual CRS GeometryCRS();

private:
	bool OpenKMLFile(std::string const & fileName);
	void CloseKMLFile();
	bool ExtractPaths();

	std::unique_ptr<char> AdvanceToNextTag();
	KMLElement GetCurrentElementValue();
	bool CompareTags(char const * tag1, char const * tag2) const;
	void ExtractNameFromKMLElement(KMLElement * element, int pathID);
	bool ExtractCoordinatesFromKMLElement(KMLElement * element, int pathID);

public:
	const FileFormat parserSupportedFormat = FileFormat::kml;

private:
	std::fstream kmlFile;

	std::string * pathsNames;

	CRS geometryCRS = CRS::WGS84; //KML is always WGS84
	unsigned int zone; //For use with UTM CRS only
	bool isNorthernHemisphere; //For use with UTM CRS only
};