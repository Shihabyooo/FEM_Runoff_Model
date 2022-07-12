//Current known limitation:
	//Assumes (and checks for, aborts otherwise) an SHP with only Polyline geometry.
	//Assumes (and checks for, warns but doens't abort) geometries with single parts (but any number of vertices, up to limit that causes length bug above).
	//Doesn't check the SHP's CRS.

//Things to add:
	//Implement and ByteToDouble() and switch uses of casting char* to double* to it instead.
	//Rename "byte" to "bytes"
	//Research and -accordingly- adjust handling of SHP files with geometries with more than 1 part.
	//Allow the use of SHP files with point geometries (Will output a single path, with vertices going from first point to last inorder).
	//Allow use of multiple geometries SHP files, in which only polylines will be considered.
	//Parse DBF files and look for names to label the paths with (stored as std::string* or char**).
	//Parse PRJ files and and save the CRS of the SHP (enum?) for use by ProfileMaker when interpolating elevations (i.e. to determine whether conversion between UTM/GCS is needed).

#pragma once
#include <fstream>
#include <iostream>
#include "Vector_File_Parser.h"

class SHPParser : public FileParser
{
public:
	SHPParser();
	~SHPParser();

	bool LoadGeometry(std::string const & fileName);
	void UnLoadGeometry();

	virtual CRS GeometryCRS();
	virtual unsigned int UTMZone();
	virtual bool IsNorthernHemisphere();

private:
	const std::string RemoveFileExtension(const std::string &fileName) const;
	bool CheckFileExistance(const std::string &fileNamePrefix) const;
	void AllocateVertsArray();
	bool LoadSHPParameters(const std::string &fileNamePrefix);	//this method will load the SHX file, and fill out pathsCount and vertsCount. Future implementation: parse DBF file to look for path names.
	bool LoadCRSDetails(const std::string &fileNamePrefix);
	bool ExtractPaths(const std::string &fileNamePrefix);
	long int BytesToInt32(const char bytes[4], bool isBigEndian) const;
	double BytesToDouble(char byte[8], bool isBigEndian) const;	//for Future implementation. For now, I'm sticking to static_cast<double> when reading the file.

public: 
	const FileFormat parserSupportedFormat = FileFormat::shapeFile;

private:
	

	std::string * pathsNames;
	long int * vertsCount;

	CRS geometryCRS = CRS::undefined;
	unsigned int zone; //For use with UTM CRS only
	bool isNorthernHemisphere; //For use with UTM CRS only
};