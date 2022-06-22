#pragma once
#include <iostream>
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
#define COLOUR_BLUE Colour(0.0f, 0.0f, 1.0f)
#define COLOUR_RED Colour(1.0f, 0.0f, 0.0f)
#define COLOUR_GREEN Colour(0.0f, 1.0f, 0.0f)
#define COLOUR_BLACK Colour()
#define COLOUR_WHITE Colour(1.0f)
#define COLOUR_GRAY Colour(0.5f)
#define COLOUR_MAGENTA Colour(1.0f, 0.0f, 1.0f)
#define COLOUR_CYAN Colour(0.0f, 1.0f, 1.0f)
#define COLOUR_YELLOW Colour(1.0f, 1.0f, 0.0f)
#define COLOUR_ORANGE Colour(1.0f, 0.5f, 0.0f)
#define COLOUR_LIME Colour(0.5f, 1.0f, 0.0f)
#define COLOUR_PURPLE Colour(0.5f, 0.0f, 1.0f)
#pragma endregion

//Mapping related defines.
#define UTM_FALSE_EASTING (double)(500000.0F)
#define UTM_FALSE_NORTHING (double)(10000000.0F)
#define UTM_MERIDIAN_SCALE (double)(0.9996F)
#define PI_CONSTANT (double)(3.14159265359F)
#define WGS84_EARTH_RADIUS_EQUATOR (double)(6378137.0F)
#define WGS_EARTH_RADIUS_POLES (double)(6356752.3142F)
#define WGS84_ELIPSOID_FLATTENING (double)(1.0F/298.257223563F)

enum class CRS
{
	WGS84,	//Geographic CRS, EPSG 4326
	UTM,	//Universal Tranverse Mercator, projected CRS
	undefined
};

enum class Solver
{
	Simple, //Compute invert of paramters matrix and multiply with the RHS vector
	GaussJordan, //Gauss-Jordan Elimination and backwards substitution
	Jacobi,
	SOR, //Successive Overrelaxation (Or Gauss-Seidel when waight = 1)
	PCG, //Preconditioned Congjugate Gradient
	BiCG, //Biconjugate Gradient
	GMRES //Generalized Minimum Residual
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
	Vector2(float const values[2]);
	
	Vector2 operator+(Vector2 const & vec2) const;
	Vector2 operator-(Vector2 const & vec2) const;
	Vector2 operator*(float const & scalar) const;
	Vector2 & operator= (Vector2 const & vec2);

	Vector2 Normalize(Vector2 const & min, Vector2 const & max) const; //Returns normalized axes relative to min and max, result is between 0.0f to 1.0f if point is inside range min-max.

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

	Vector2 Normalize(Vector2Int const & min, Vector2Int const & max) const; //Returns normalized axes relative to min and max, result is between 0.0f to 1.0f if point is inside range min-max.

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

enum class LogEntryType
{
	normal, warning, error, success
};

struct LogEntry
{
public:
	LogEntry(std::string const & _content, LogEntryType _type)
	{
		content = _content;
		type = _type;
	};

	//TODO add timestamp
	std::string content;
	LogEntryType type;
};

struct Colour
{
public:
	Colour()
	{
		r = g = b = 0.0f;
		a = 1.0f;
	}

	Colour(float grayShade)
	{
		r = g = b = grayShade;
		a = 1.0f;
	}

	Colour(float _r, float _g, float _b)
	{
		r = _r;
		g = _g;
		b = _b;
		a = 1.0f;
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

//Helper functions
//TODO research whether there is any difference between having them on source file vs header file.

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

static inline void Print(Vector2 const & vec)
{
	std::cout << vec.x << ", " << vec.y << "\n";
}

static inline void Print(Vector2D const & vec)
{
	std::cout << vec.x << ", " << vec.y << "\n";
}

//Projection, Coordinate transformations and interpolation
std::unique_ptr<double> ToUTM(double lng, double lat);
std::unique_ptr<double> ToWGS84(double easting, double northing, bool isNortherHemisphere, int zone);

Vector2D ProjectPoint(Vector2D const & point); //converts from geographic to UTM.point.x = longitude, point.y = latitude.\
												result.x = easting, result.y = northing.
Vector2D WrapPoint(Vector2D const & point, bool isNorthernHemisphere, int zone); //Converts from UTM to geographic. point.x = easting, point.y = northing.\
																				 result.x = longitude, result.y = latitude


double LinearInterpolationNormalized(double normalizedPoint, double values[2]);
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