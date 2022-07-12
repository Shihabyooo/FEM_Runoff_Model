//This, plus SHP_Parser and KML_Parser were originally part of another project, copied here with minor changes.

#pragma once
#include "Globals.hpp"
//#include "MatricesPP.hpp"

enum class FileFormat
{
	shapeFile, kml, csv, unsupported
};

class FileParser
{
public:
	FileParser() {};
	~FileParser() {};

	virtual bool LoadGeometry(std::string const & fileName) { return false; };
	virtual void UnLoadGeometry() {};

	virtual Matrix_f64 const * const GetPathByID(int id)
	{
		if (id >= pathsCount || !isPathLoaded)
			return NULL;
		else
			return &verts[id];
	};

	virtual bool IsPathLoaded() { return isPathLoaded; };

	virtual CRS GeometryCRS() { return geometryCRS; };
	virtual unsigned int UTMZone() { return 0; };
	virtual bool IsNorthernHemisphere() { return true; };

public:
	const FileFormat parserSupportedFormat = FileFormat::unsupported;

protected:
	bool isPathLoaded = false;
	long int pathsCount = 0;
	Matrix_f64 * verts;
	CRS geometryCRS = CRS::undefined;
	unsigned int zone; //For use with UTM CRS only
	bool isNorthernHemisphere; //For use with UTM CRS only
};