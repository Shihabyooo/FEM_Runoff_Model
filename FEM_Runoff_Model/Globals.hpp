#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <MatricesPP.hpp>

//TODO alot of the variables, defs, structs and methods of this module were implemented for certain features that are no longer\
part of this program or were dropped. A lot of cleanup is required here.

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

enum class SpatialSamplingMethod
{
	average = 0,
	median = 1,
	majority = 2,
	nearest = 3,
};

enum class TimeUnit
{
	second = 0,
	minute = 1,
	hour = 2,
	day = 3
};

enum class LossModel
{
	none = 0,
	initialConst = 1,
	scsCN = 2,
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
	bool operator== (Vector2D const & vec2) const;

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
	bool ContainsInclusive(Vector2D point); //includes edges

	Vector2D minCorner = Vector2D();
	Vector2D maxCorner = Vector2D();
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

	Rect xy = Rect();
	double z_SW = 0.0, z_SE = 0.0, z_NE = 0.0, z_NW = 0.0;
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
	void ComputeSum();
	double HoursToLocalUnits(double time) const;
	double SampleRate(double timeSinceStart) const; //timeSinceStart in hours. Returns rate as mm/hr
	double SampleCummulativePreciptation(double timeSinceStart) const; //timeSinceStart in hours. Returns cummulative (total) precipitation up to that time in mm.

	size_t size = 0;
	TimeUnit timeUnit = TimeUnit::hour;
	double precipitationSum = 0.0; 
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
	//InterpolationType topographySamplingMethod = InterpolationType::nearest;

	//TODO expose to GUI
	size_t outletNode = 0; //the watershed's outlet. 
	//Precipitation
	bool variablePrecipitation = true; //false use unitTimeseries for all elements for all periods.
										//true, use gridded precipitation of 1D time series.

	//bool griddedPrecipitation = false; //false: precipitation value from 1D time series for all elements. unitTimeSeries must be set to TS (must be pair<double, double>[]).
									//True: precipitationRastersFolder must be set
	
	//std::string precipitationRastersFolder = ""; //Folder must contain rasters named appropriatly for each time step of simulation. To be implemented.

	TimeSeries unitTimeSeries;

	//Temporal intpolertion is unused in current implmenetation.
	//InterpolationType precipitationTemporalInterpolationType = InterpolationType::linear; //should either be linear or cubic. Nearest should never 
																				//be used unless timeSeries resolution is close to or finer than
																				//simulation timeStep
	//InterpolationType precipitationSpatialInterpolationType = InterpolationType::nearest;

	LossModel lossModel = LossModel::none;

	void * lossModelParams = NULL;
	unsigned int scsCN = 0; //This is placeholder, untill proper solution that supports gridded SCS is implemented.

	//Hydraulic Parameters
	bool variableManningCoefficients = false; //If false: using fixedManningCoeffients for all elements. If true: gridded manningCoefficientsRaster must be supplied.
	double fixedManningCoeffient = -1.0f; //must be positive value > 0.0
	std::string manningCoefficientRasterPath = "";

	//temporal params
	double timeStep = 0.5; //delta T, in hours. e.g. 0.5 = 30 minutes, 1.0 = 1 hour.
	double startTime = 0.0; //should be left at 0.0
	double endTime = 10.0f; //hours after startTime to end simulation.

	//FEM related params
	bool useLumpedForm = true; //if false, uses consistant formulation
	//femOmega is a weighting factor to control temporal approximation.
	//0.0 = Forward difference, 1.0 = Backward difference 0.5 = Central difference (Crank-Nicholson method)
	double femOmega = 0.5; 

	//Solver related params
	Solver solverType = Solver::Auto;
	double residualThreshold = -1.0; //Negative value -> use default threshold. Only for iterative solvers.
	double weight = -1.0; //Negative value -> use default weight. Only for weighted solvers.
	size_t maxIterations = 0; //0 -> Use default value. Only for iterative solvers.
	double externalResidualTreshold = 0.00001;
	size_t maxExternalIterations = 1000; //for internal loop.
};

struct MeshGeneratorParameters
{
public:
	std::vector<Vector2D> const * boundary = NULL;
	size_t resolution = 10; //must be greater than 2.
	double internalPadding = 0.001; //must be real, positive value.
	double rayCastPadding = 5.0; //must be real, positive value.
};

struct InitialAndConstantParams
{
public:
	double initialLoss = 0.0; //in mm
	double constRate = 0.0; //in mm/hr
};

//Helper functions
template<typename T>
static T Min(T const & a, T const & b)
{
	return (a > b ? b : a);
}

template<typename T>
static T Max(T const & a, T const & b)
{
	return (a > b ? a : b);
}

template<typename T>
static T Clamp(T const & a, T const & b, T const & c) //Clamps a to range b, c.
{
	T min = Min(b, c);
	T max = Max(b, c);
	return (a > max ? max : (a < min ? min : a));
}

template<typename T>
static T Average(std::vector<T> const & values)
{
	T result = 0.0;
	
	for (auto it = values.begin(); it != values.end(); ++it)
		result += *it;
	
	return result / static_cast<T>(values.size());
}

template<typename T>
static T Median(std::vector<T> const & values)
{
	if (values.size() < 1)
		return T();
	else if (values.size() < 2)
		return values[0];

	std::vector<T> orderedValues = values;
	std::sort(orderedValues.begin(), orderedValues.end());

	if (orderedValues.size() % 2 == 1)
		return orderedValues[(orderedValues.size() - 1) / 2];
	
	T avg = orderedValues[orderedValues.size() / 2] + orderedValues[(orderedValues.size() / 2) - 1];
	return avg / 2.0; //This may cause issues...
}

template<typename T>
static T Majority(std::vector<T> const & values) //Returns smallest element if no majority was found. Returns majority with smallest value for T if multiple had equal frequency.
{
	if (values.size() < 1)
		return T();
	else if (values.size() < 2)
		return values[0];

	std::vector<T> orderedValues = values;
	std::sort(orderedValues.begin(), orderedValues.end());

	std::vector<std::pair<T, size_t>> frequency;
	frequency.push_back(std::pair < T, size_t>(orderedValues[0], 1));

	for (size_t i = 1; i < values.size(); i++)
	{
		if (orderedValues[i] == frequency.back().first)
			frequency.back().second++;
		else
			frequency.push_back(std::pair < T, size_t>(orderedValues[i], 1));
	}

	T result;
	size_t maxFrequency = 0;
	for (auto it = frequency.begin(); it != frequency.end(); ++it)
		if (it->second > maxFrequency)
		{
			result = it->first;
			maxFrequency = it->second;
		}

	return result;
}

template<typename T>
static T Majority(std::vector<T> const & values, double tolerance) //for floating points. Tolerence is the range within two values are assumed equal.
{
	std::vector<T> adjustedValues = values;

	if (tolerance > 0.0)
	{
		std::sort(adjustedValues.begin(), adjustedValues.end());

		for (size_t i = 1; i < adjustedValues.size(); i++)
		{
			if (adjustedValues[i] - adjustedValues[i - 1] < tolerance)
				adjustedValues[i] = adjustedValues[i - 1];
		}
	}

	return Majority(adjustedValues);
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
