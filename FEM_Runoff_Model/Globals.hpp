#pragma once
#include <iostream>
#include <vector>
#include <MatricesPP.hpp>

#define PROGRAM_NAME "FEM_Runoff_Model"

//Error codes.
#define SUCCESS 0
#define FAILED_MAIN_WINDOW_INITIALIZE 1
#define FAILED_VIEWPORT_CREATE 2

//Log Types definitions.
#define LOG_SUCCESS	LogEntryType::success
#define LOG_ERROR	LogEntryType::error
#define LOG_NORM	LogEntryType::normal
#define LOG_WARN	LogEntryType::warning

//Colours
#pragma region Colour templates
#define COLOUR_BLUE Colour(0.0, 0.0, 1.0)
#define COLOUR_RED Colour(1.0, 0.0, 0.0)
#define COLOUR_GREEN Colour(0.0, 1.0, 0.0)
#define COLOUR_BLACK Colour()
#define COLOUR_WHITE Colour(1.0)
#define COLOUR_GRAY Colour(0.5)
#define COLOUR_MAGENTA Colour(1.0, 0.0, 1.0)
#define COLOUR_CYAN Colour(0.0, 1.0, 1.0)
#define COLOUR_YELLOW Colour(1.0, 1.0, 0.0)
#define COLOUR_ORANGE Colour(1.0, 0.5, 0.0)
#define COLOUR_LIME Colour(0.5, 1.0, 0.0)
#define COLOUR_PURPLE Colour(0.5, 0.0, 1.0)
#pragma endregion

//Mapping related defines.
#define UTM_FALSE_EASTING (double)(500000.0)
#define UTM_FALSE_NORTHING (double)(10000000.0)
#define UTM_MERIDIAN_SCALE (double)(0.9996)
#define PI_CONSTANT (double)(3.14159265359)
#define WGS84_EARTH_RADIUS_EQUATOR (double)(6378137.0)
#define WGS_EARTH_RADIUS_POLES (double)(6356752.3142)
#define WGS84_ELIPSOID_FLATTENING (double)(1.0F/298.257223563)

//Misc Defines
//#define OUT_FILE_PATH "output.csv"
#define LOG_FILE_PATH "log_file.txt" //in same dir as exe
#define LOG_TO_CLI //also outputs the logged messages to the CLI (using std::cout)

enum class CRS
{
	WGS84,	//Geographic CRS, EPSG 4326
	UTM,	//Universal Tranverse Mercator, projected CRS
	undefined
};

enum class Solver
{
	Auto = 0,
	Simple = 1, //Compute invert of paramters matrix and multiply with the RHS vector
	Gaussian = 2, //Gauss Elimination and backwards substitution
	Jacobi = 3, //(Weighted) Jacobi
	SOR = 4, //Successive Overrelaxation (Or Gauss-Seidel when waight = 1)
	PCG = 5, //Preconditioned Congjugate Gradient
	BiCG = 6, //(Preconditioned) Biconjugate Gradient
	CGS = 7, //(Preconditioned) Conjugate Gradient Squared.
	//GMRES = 8 //Generalized Minimal Residual
};

enum class LogEntryType
{
	normal, warning, error, success
};

enum class InterpolationType
{
	nearest = 0,
	linear = 1, //linear for 1D, bilinear for 2D
	cubic = 2 //cubic for 1D, bicbic for 2D
};

enum class TimeUnit
{
	second = 0,
	minute = 1,
	hour = 2,
	day = 3
};

enum class ElementType
{	
	triangle = 0,
	rectangle = 1,
	undefined = 2
};

struct Vector2;
struct Vector2Int;
struct Vector2D;

struct Vector2
{
public:
	Vector2();
	Vector2(float _x, float _y);
	Vector2(Vector2 const & vec);
	Vector2(Vector2Int const & vec);
	explicit Vector2(Vector2D const & vec);
	Vector2(float const values[2]);
	
	Vector2 operator+(Vector2 const & vec2) const;
	Vector2 operator-(Vector2 const & vec2) const;
	Vector2 operator*(float const & scalar) const;
	Vector2 & operator= (Vector2 const & vec2);

	Vector2 Normalize(Vector2 const & min, Vector2 const & max) const; //Returns normalized axes relative to min and max, result is between 0.0 to 1.0 if point is inside range min-max.
	float DistanceTo(Vector2 const & vec2) const;
	bool WithinCircle(Vector2 const & centre, float radius) const;

	float x;
	float y;
};

struct Vector2Int
{
public:
	Vector2Int();
	Vector2Int(int _x, int _y);
	Vector2Int(int const values[2]);

	Vector2Int operator+(Vector2Int const & vec2) const;
	Vector2Int operator-(Vector2Int const & vec2) const;
	Vector2Int operator*(int const & scalar) const;
	Vector2Int & operator= (Vector2Int const & vec2);
	bool operator== (Vector2Int const & vec2);

	Vector2 Normalize(Vector2Int const & min, Vector2Int const & max) const; //Returns normalized axes relative to min and max, result is between 0.0 to 1.0 if point is inside range min-max.

	int x;
	int y;
};

struct Vector2D
{
public:
	Vector2D();
	Vector2D(double _x, double _y);
	Vector2D(Vector2D const & vec);
	Vector2D(Vector2 const & vec);
	Vector2D(Vector2Int const & vec);
	Vector2D(double const values[2]);

	Vector2D operator+(Vector2D const & vec2) const;
	Vector2D operator-(Vector2D const & vec2) const;
	Vector2D operator*(double const & scalar) const;
	Vector2D & operator= (Vector2D const & vec2);

	Vector2D Normalize(Vector2D const & min, Vector2D const & max) const; //Returns normalized axes relative to min and max, result is between 0.0F to 1.0F if point is inside range min-max.
	double DistanceTo(Vector2D const & vec2) const;
	bool WithinCircle(Vector2D const & centre, double radius) const;

	double x;
	double y;
};

struct Rect
{
public:
	Rect();
	Rect(Rect const & rect);
	Rect(Vector2D const & cornerSW, Vector2D const & cornerNE);
	Rect(Vector2D const & cornerSW, double width, double height);
	Rect(double width, double height, Vector2D const & cornerNE);

	Rect & operator= (Rect const & rect);

	double Width() const;
	double Height() const;
	double Area() const;
	void Translate(Vector2D delta);
	void Translate(Vector2D direction, double magnitude);
	bool Contains(Vector2D point); //not including edges
	
	Vector2D minCorner;
	Vector2D maxCorner;
};

struct Rect3D
{
public:
	Rect3D();
	Rect3D(Rect3D const & rect3D);
	Rect3D(Rect const & xyRect);
	Rect3D(Rect const & xyRect, double uniformHeight);
	Rect3D(Rect const & xyRect, double heightSW, double heightSE, double heightNE, double heightNW);
	Rect3D(Rect const & xyRect, double const heights[4]); //heights: SW->SE->NW->NE
	
	Rect3D & operator= (Rect3D const & rect3D);

	bool Contains(Vector2D const & planarPos);
	Vector2D MinCorner() const;
	Vector2D MaxCorner() const;

	Rect xy;
	double z_SW, z_SE, z_NE, z_NW;
};

//Name should probably be Grid3x3. When I thought of this struct I had vertices as the defining point (pun not intended) rather tha cells.
struct Grid4x4 //Grid of 4 by 4 cells, SW-anchored, each cell with same dimension (but height can be different than width)
{
public:
	Grid4x4();
	Grid4x4(Grid4x4 const & grid);
	Grid4x4(Rect const & centralRect);
	Grid4x4(Rect const & centralRect, double const uniformZValue);
	Grid4x4(Rect const & centralRect, double const zValues[16]); //zValue: From east to west, south to north
	Grid4x4(Rect const & centralRect, double const zValues[4][4]); //zValue = double[4][4], from east to west, south to north
	Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight);
	Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight, double const uniformZValue);
	Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight, double const zValues[16]); //zValue = double[16], from east to west, south to north
	Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight, double const zValues[4][4]); //zValue = double[4][4], from east to west, south to north

	Grid4x4 & operator= (Grid4x4 const & grid);

	void GridFromCentralRect(Rect const & centralRect);
	void GridFromSWCorner(Vector2D const & cornerSW, double const cellWidth, double const cellHeight);
	void SetUniformZValue(double const zValue);
	void SetZValues(double const zValues[16]); //zValue = double[16], from west to east, south to north
	void SetZValues(double const zValues[4][4]); //zValue = double[4][4], from west to east, south to north

	double CellWidth() const;
	double CellHeight() const;

	Vector2D GridNode(unsigned int i, unsigned int j) const; //i, j < 4, otherwise returns a zeroed Vector2D.
	Rect3D GetCell(unsigned int i, unsigned int j) const; //i, j < 3. Anchor is SW cell, order west to east, south to north. OOB Returns zeroed Rect3D
	
	double x[4];
	double y[4];
	double z[4][4];
};

struct LogEntry
{
public:
	LogEntry(std::string const & _content, LogEntryType _type)
	{
		content = _content;
		type = _type;
	};

	std::string content; //content is prefexed with a timestamp.
	LogEntryType type;
};

struct Colour
{
public:
	Colour()
	{
		r = g = b = 0.0;
		a = 1.0;
	}

	Colour(float grayShade)
	{
		r = g = b = grayShade;
		a = 1.0;
	}

	Colour(float _r, float _g, float _b)
	{
		r = _r;
		g = _g;
		b = _b;
		a = 1.0;
	}

	Colour(float _r, float _g, float _b, float _a)
	{
		r = _r;
		g = _g;
		b = _b;
		a = _a;
	}

	float const * Array() const
	{
		float array[4]{ r, g, b, a };
		return array;
	}

	float r, g, b, a;
};

struct TimeSeries //Always ensures at least a time series of 2 entries.
{
public:
	TimeSeries();
	TimeSeries(size_t _size); //Will force size to 2 if _size < 2
	TimeSeries(std::vector<std::pair<size_t, double>> const & ts);
	~TimeSeries();

	void operator= (TimeSeries const & ts);

	bool IsValid() const;
	void AdjustSize(size_t newSize);
	double HoursToLocalUnits(double time) const;
	double SampleRate(double timeSinceStart) const;// , double timeSpan, InterpolationType interpolationType) const; //timeSinceStart in hours. Returns rate as mm/hr

	size_t size = 0;
	TimeUnit timeUnit = TimeUnit::hour;
	//probably would make more sense to use a vector...
	std::pair<size_t, double> * series = NULL;
};

struct ModelParameters
{
	//NOTE! Nodes and Elements should be already loaded and set by the time this struct is needed, so they are not included here, but still
	//critical for model.
public:
	ModelParameters();	
	~ModelParameters();

	std::string modelName = "untitled_model";

	//Topographic/geogrpahic params
	std::string demPath = "";
	std::string slopesPath = "";
	std::string fdrPath = "";
	InterpolationType topographySamplingMethod = InterpolationType::nearest;

	//TODO expose to GUI
	size_t outletNode = 0; //the watershed's outlet. 
	//Precipitation
	bool variablePrecipitation = true; //false use unitTimeseries for all elements for all periods.
										//true, use gridded precipitation of 1D time series.

	//bool griddedPrecipitation = false; //false: precipitation value from 1D time series for all elements. unitTimeSeries must be set to TS (must be pair<double, double>[]).
									//True: precipitationRastersFolder must be set
	
	//std::string precipitationRastersFolder = ""; //Folder must contain rasters named appropriatly for each time step of simulation. To be implemented.

	TimeSeries unitTimeSeries;	//(Discription bellow is from old impl that was a simple std::pair<double, double> array
								//first double is time relative to startTime, second is in incremental mm.
								//e.g. if simulation startTime 0.0 is equivalent to Jan 1st 00:00 AM, the TS:
								//<0.0, 0.0>	-> 0mm at start. Should always be the case
								//<1.0, 5.0>	->  5mm between 00:00 AM and 01:00 AM
								//<2.0, 15.0>	-> 15mm between 01:00 AM and 02:00 AM
								//<3.0, 30.0>	-> 30mm between 02:00 AM and 03:00 AM
								//<4.0, 10.0>	-> 10mm between 03:00 AM and 04:00 AM	
								//<10.0, 7.5>	-> 7.5mm between 04:00 AM and 010:00 AM	
								//Any precipitation after 10AM will be assumed zero.
								//Preciptation between 

	InterpolationType precipitationTemporalInterpolationType = InterpolationType::linear; //should either be linear or cubic. Nearest should never 
																				//be used unless timeSeries resolution is close to or finer than
																				//simulation timeStep
	InterpolationType precipitationSpatialInterpolationType = InterpolationType::nearest;

	//double fixedPrecipitationValue = -1.0f; //must be positive value > 0.0

	//Hydraulic Parameters
	bool variableManningCoefficients = false; //If false: using fixedManningCoeffients for all elements.
												//If true: gridded manningCoefficientsRaster must be supplied.

	double fixedManningCoeffient = -1.0f; //must be positive value > 0.0
	std::string manningCoefficientRasterPath = "";

	//bool useBuiltInLossModel = true;
	//bool useHydrologicClassGrid = false; //if true, hydrologic class raster must be set
	//std::string hydrologicClassRaster;

	//temporal params
	double timeStep = 0.5; //delta T, in hours. e.g. 0.5 = 30 minutes, 1.0 = 1 hour.
	double startTime = 0.0; //should be left at 0.0
	double endTime = 10.0f; //hours after startTime to end simulation.

	//FEM related params
	bool useLumpedForm = true; //if false, uses consistant formulation
	double femOmega = 0.5; //A weighting factor to control temporal approximation. 0.0 = Forward difference, 1.0 = Backward difference\
							0.5 = Central difference (Crank-Nicholson method)

	//Solver related params
	//ElementType meshType = ElementType::undefined;
	Solver solverType = Solver::Auto;
	double residualThreshold = -1.0; //Negative value -> use default threshold. Only for iterative solvers.
	double weight = -1.0; //Negative value -> use default weight. Only for weighted solvers.
	size_t maxIterations = 0; //0 -> Use default value. Only for iterative solvers.
	double internalResidualTreshold = 0.00001;
	size_t maxInternalIterations = 1000; //for internal loop.
	
	//output
	int nthDurationToOutput = 1; //as in output every nth simulation frame to disk. >=1 means output all simulated from. 2 means output every other frame, etc.
	//TODO implement output formating.
};

struct MeshGeneratorParameters
{
public:
	ElementType meshType = ElementType::undefined; //Required Triangle or Rectangle
	std::vector<Vector2D> const * boundary = NULL; //Required except if meshType == Triangle && useCustomNodes == true.

	//for triangular elements
	bool useCustomNodes = true; //If true, nodesList must be supplied (to user, inNodesListPath is what must be supplied).
	std::vector<Vector2D> * inNodesList = NULL;
	std::string inNodesListPath = "";
	Vector2D * outSuperTriangleNodes = NULL; //optional. If supplied, must point to an array of 6 Vector2Ds.
	double superTrianglePadding = -1.0; // must be positive real number greater than 0.0

	//for rectangular elements
	size_t resolution = 10; //must be greater than 2.
	double internalPadding = 0.001; //must be real, positive value.
	double rayCastPadding = 1.0; //must be real, positive value.
};

//Helper functions

static double Min(double const & a, double const & b)
{
	return (a > b ? b : a);
}

static float Min(float const & a, float const & b)
{
	return (a > b ? b : a);
}

static int Min(int const & a, int const & b)
{
	return (a > b ? b : a);
}

static size_t Min(size_t const & a, size_t const & b)
{
	return (a > b ? b : a);
}

static double Max(double const & a, double const & b)
{
	return (a > b ? a : b);
}

static float Max(float const & a, float const & b)
{
	return (a > b ? a : b);
}

static int Max(int const & a, int const & b)
{
	return (a > b ? a : b);
}

static size_t Max(size_t const & a, size_t const & b)
{
	return (a > b ? a : b);
}

static float Clamp(float const & a, float const & b, float const & c)
{
	float min = Min(b, c);
	float max = Max(b, c);
	return (a > max ? max : (a < min ? min : a));
}

static int Clamp(int const & a, int const & b, int const & c)
{
	int min = Min(b, c);
	int max = Max(b, c);
	return (a > max ? max : (a < min ? min : a));
}

static inline void Print(Vector2 const & vec, bool sameLine = false)
{
	std::cout << vec.x << ", " << vec.y;
	
	if (sameLine)
		std::cout << " ";
	else
		std::cout << "\n";
}

static inline void Print(Vector2D const & vec, bool sameLine = false)
{
	std::cout << vec.x << ", " << vec.y;

	if (sameLine)
		std::cout << " ";
	else
		std::cout << "\n";
}

//Projection, Coordinate transformations and interpolation
std::unique_ptr<double> ToUTM(double lng, double lat);
std::unique_ptr<double> ToWGS84(double easting, double northing, bool isNortherHemisphere, int zone);

Vector2D ProjectPoint(Vector2D const & point); //converts from geographic to UTM.point.x = longitude, point.y = latitude.\
												result.x = easting, result.y = northing.
Vector2D WrapPoint(Vector2D const & point, bool isNorthernHemisphere, int zone); //Converts from UTM to geographic. point.x = easting, point.y = northing.\
																				 result.x = longitude, result.y = latitude


double LinearInterpolationNormalized(double normalizedPoint, double const values[2]);
double CubicInterpolationNormalized(double const normalizedPoint, double const values[4]); 

double BilinearInterpolation(Vector2D const & point, Grid4x4 const & grid);
double BilinearInterpolation(Vector2D const & point, Rect3D const & grid);
double BicubicInterpolation(Vector2D const & point, Grid4x4 const & grid);

//double BicubicInterpolation(double const x, double const y,
//							double const x0, double const y0, double const z0,
//							double const x1, double const y1, double const z1,
//							double const x2, double const y2, double const z2,
//							double const x3, double const y3, double const z3);
//
//double BicubicInterpolation(double const x, double const y, double const * boundsX, double const * boundsY, double const ** boundsZ);