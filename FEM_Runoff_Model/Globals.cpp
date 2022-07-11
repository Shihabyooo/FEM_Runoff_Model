#include "Globals.hpp"

#pragma region Vector2 Defintions

Vector2::Vector2()
{
	x = 0.0f;
	y = 0.0f;
}

Vector2::Vector2(float _x, float _y)
{
	x = _x;
	y = _y;
}

Vector2::Vector2(Vector2 const & vec)
{
	x = vec.x;
	y = vec.y;
}

Vector2::Vector2(Vector2Int const & vec)
{
	x = static_cast<float>(vec.x);
	y = static_cast<float>(vec.y);
}

Vector2::Vector2(Vector2D const & vec)
{
	x = static_cast<float>(vec.x);
	y = static_cast<float>(vec.y);
}

Vector2::Vector2(float const values[2])
{
	x = values[0];
	y = values[1];
}

Vector2 Vector2::operator+(Vector2 const & vec2) const
{
	return (Vector2(x + vec2.x, y + vec2.y));
}

Vector2 Vector2::operator-(Vector2 const & vec2) const
{
	return (Vector2(x - vec2.x, y - vec2.y));
};

Vector2 Vector2::operator*(float const & scalar) const
{
	return (Vector2(x * scalar, y * scalar));
};

Vector2 & Vector2::operator= (Vector2 const & vec2)
{
	x = vec2.x;
	y = vec2.y;
	return *this;
}

Vector2 Vector2::Normalize(Vector2 const & min, Vector2 const & max) const
{
	Vector2 delta = max - min;
	return Vector2(	(x - min.x) / delta.x,
					(y - min.y) / delta.y);
}

float Vector2::DistanceTo(Vector2 const & vec2) const
{
	return sqrt(pow(x - vec2.x, 2) + pow(y - vec2.y, 2));
}

bool Vector2::WithinCircle(Vector2 const & centre, float radius) const
{
	return ((x - centre.x) * (x - centre.x) + (y - centre.y) * (y - centre.y)) < (radius * radius);
}

#pragma endregion

#pragma region Vector2 Defintions

Vector2Int::Vector2Int()
{
	x = 0;
	y = 0;
}

Vector2Int::Vector2Int(int _x, int _y)
{
	x = _x;
	y = _y;
}

Vector2Int::Vector2Int(int const values[2])
{
	x = values[0];
	y = values[1];
}

Vector2Int Vector2Int::operator+(Vector2Int const & vec2) const
{
	return (Vector2Int(x + vec2.x, y + vec2.y));
}

Vector2Int Vector2Int::operator-(Vector2Int const & vec2) const
{
	return (Vector2Int(x - vec2.x, y - vec2.y));
};

Vector2Int Vector2Int::operator*(int const & scalar) const
{
	return (Vector2Int(x * scalar, y * scalar));
};

Vector2Int & Vector2Int::operator= (Vector2Int const & vec2)
{
	x = vec2.x;
	y = vec2.y;
	return *this;
}

bool Vector2Int::operator==(Vector2Int const & vec2)
{
	return (x == vec2.x) && (y == vec2.y);
}

Vector2 Vector2Int::Normalize(Vector2Int const & min, Vector2Int const & max) const
{
	Vector2 delta = max - min;
	return Vector2(	(x - min.x) / delta.x,
					(y - min.y) / delta.y);
}

#pragma endregion

#pragma region Vector2D Definitions

Vector2D::Vector2D()
{
	x = 0.0;
	y = 0.0;
}

Vector2D::Vector2D(double _x, double _y)
{
	x = _x;
	y = _y;
}

Vector2D::Vector2D(Vector2D const & vec)
{
	x = vec.x;
	y = vec.y;
}

Vector2D::Vector2D(Vector2 const & vec)
{
	x = static_cast<double>(vec.x); //unncessary? Compiler alread handles casting well from smaller types to doubles.
	y = static_cast<double>(vec.y);
}

Vector2D::Vector2D(Vector2Int const & vec)
{
	x = static_cast<double>(vec.x);
	y = static_cast<double>(vec.y);
}

Vector2D::Vector2D(double const values[2])
{
	x = values[0];
	y = values[1];
}

Vector2D Vector2D::operator+(Vector2D const & vec2) const
{
	return (Vector2D(x + vec2.x, y + vec2.y));
}

Vector2D Vector2D::operator-(Vector2D const & vec2) const
{
	return (Vector2D(x - vec2.x, y - vec2.y));
};

Vector2D Vector2D::operator*(double const & scalar) const
{
	return (Vector2D(x * scalar, y * scalar));
};

Vector2D & Vector2D::operator= (Vector2D const & vec2)
{
	x = vec2.x;
	y = vec2.y;
	return *this;
}

Vector2D Vector2D::Normalize(Vector2D const & min, Vector2D const & max) const
{
	Vector2D delta = max - min;
	return Vector2D((x - min.x) / delta.x,
					(y - min.y) / delta.y);
}

double Vector2D::DistanceTo(Vector2D const & vec2) const
{
	return sqrt(pow(x - vec2.x , 2) + pow( y - vec2.y , 2));
}

bool Vector2D::WithinCircle(Vector2D const & centre, double radius) const
{
	return ((x - centre.x) * (x - centre.x) + (y - centre.y) * (y - centre.y)) < (radius * radius);
}

#pragma endregion

#pragma region Rect Definitions

Rect::Rect()
{
	minCorner = Vector2D();
	maxCorner = Vector2D();
}

Rect::Rect(Rect const & rect)
{
	*this = rect;
}

Rect::Rect(Vector2D const & cornerSW, Vector2D const & cornerNE)
{
	minCorner = cornerSW;
	maxCorner = cornerNE;
}

Rect::Rect(Vector2D const & cornerSW, double width, double height)
{
	minCorner = cornerSW;
	maxCorner = cornerSW + Vector2D(width, height);
}

Rect::Rect(double width, double height, Vector2D const & cornerNE)
{
	maxCorner = cornerNE;
	minCorner = cornerNE - Vector2D(width, height);
}

Rect & Rect::operator=(Rect const & rect)
{
	minCorner = rect.minCorner;
	maxCorner = rect.maxCorner;
	return *this;
}

double Rect::Width() const
{
	return maxCorner.x - minCorner.x;
}

double Rect::Height() const
{
	return maxCorner.y - minCorner.y;
}

double Rect::Area() const
{
	return Width() * Height();
}

void Rect::Translate(Vector2D delta)
{
	minCorner = minCorner + delta;
	maxCorner = maxCorner + delta;
}

void Rect::Translate(Vector2D direction, double magnitude)
{
	Vector2D delta = direction * magnitude;
	Translate(delta);
}

bool Rect::Contains(Vector2D point)
{
	return	point.x > minCorner.x && point.x < maxCorner.x &&
			point.y > minCorner.y && point.y < maxCorner.y;
			
}

#pragma endregion

#pragma region Rect3D Defitions

Rect3D::Rect3D()
{
	xy = Rect();
	z_SW = z_SE = z_NE = z_NW = 0.0;
}

Rect3D::Rect3D(Rect3D const & rect3D)
{
	*this = rect3D;
}

Rect3D::Rect3D(Rect const & xyRect)
{
	xy = xyRect;
	z_SW = z_SE = z_NE = z_NW = 0.0;
}

Rect3D::Rect3D(Rect const & xyRect, double uniformHeight)
{
	xy = xyRect;
	z_SW = z_SE = z_NE = z_NW = uniformHeight;
}

Rect3D::Rect3D(Rect const & xyRect, double heightSW, double heightSE, double heightNE, double heightNW)
{
	xy = xyRect;
	z_SW = heightSW;
	z_SE = heightSE;
	z_NE = heightNE;
	z_NW = heightNW;
}

Rect3D::Rect3D(Rect const & xyRect, double const heights[4])
{
	xy = xyRect;
	z_SW = heights[0];
	z_SE = heights[1];
	z_NW = heights[2];
	z_NE = heights[3];
}

Rect3D & Rect3D::operator=(Rect3D const & rect3D)
{
	xy = rect3D.xy;
	z_SW = rect3D.z_SW;
	z_SE = rect3D.z_SE;
	z_NE = rect3D.z_NE;
	z_NW = rect3D.z_NW;

	return *this;
}

bool Rect3D::Contains(Vector2D const & planarPos)
{
	return xy.Contains(planarPos);
}

Vector2D Rect3D::MinCorner() const
{
	return xy.minCorner;
}

Vector2D Rect3D::MaxCorner() const
{
	return xy.maxCorner;
}

#pragma endregion

#pragma region Grid4x4 Definitions

Grid4x4::Grid4x4()
{
	x[0] = x[1] = x[2] = x[3] = 0.0;
	y[0] = y[1] = y[2] = y[3] = 0.0;
	
	SetUniformZValue(0.0);
}

Grid4x4::Grid4x4(Grid4x4 const & grid)
{
	*this = grid;
}

Grid4x4::Grid4x4(Rect const & centralRect)
{
	GridFromCentralRect(centralRect);
	SetUniformZValue(0.0);
}

Grid4x4::Grid4x4(Rect const & centralRect, double const uniformZValue)
{
	GridFromCentralRect(centralRect);
	SetUniformZValue(uniformZValue);
}

Grid4x4::Grid4x4(Rect const & centralRect, double const zValues[16])
{
	GridFromCentralRect(centralRect);
	SetZValues(zValues);
}

Grid4x4::Grid4x4(Rect const & centralRect, double const zValues[4][4])
{
	GridFromCentralRect(centralRect);
	SetZValues(zValues);
}

Grid4x4::Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight)
{
	GridFromSWCorner(cornerSW, cellWidth, cellHeight);
	SetUniformZValue(0.0);
}

Grid4x4::Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight, double const uniformZValue)
{
	GridFromSWCorner(cornerSW, cellWidth, cellHeight);
	SetUniformZValue(uniformZValue);
}

Grid4x4::Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight, double const zValues[16])
{
	GridFromSWCorner(cornerSW, cellWidth, cellHeight);
	SetZValues(zValues);
}

Grid4x4::Grid4x4(Vector2D const & cornerSW, double const cellWidth, double const cellHeight, double const zValues[4][4])
{
	GridFromSWCorner(cornerSW, cellWidth, cellHeight);
	SetZValues(zValues);
}

Grid4x4 & Grid4x4::operator=(Grid4x4 const & grid)
{
	for (int i = 0; i < 4; i++)
	{
		x[i] = grid.x[i];
		y[i] = grid.y[i];

		for (int j = 0; j < 4; j++)
			z[i][j] = grid.z[i][j];
	}
}

void Grid4x4::GridFromCentralRect(Rect const & centralRect)
{
	double width = centralRect.Width();
	double height = centralRect.Height();

	x[0] = centralRect.minCorner.x - width;
	x[1] = centralRect.minCorner.x;
	x[2] = centralRect.maxCorner.x;
	x[3] = centralRect.maxCorner.x + width;

	y[0] = centralRect.minCorner.y - height;
	y[1] = centralRect.minCorner.y;
	y[2] = centralRect.maxCorner.y;
	y[3] = centralRect.maxCorner.y + height;
}

void Grid4x4::GridFromSWCorner(Vector2D const & cornerSW, double const cellWidth, double const cellHeight)
{
	for (int i = 0; i < 4; i++)
	{
		x[i] = cornerSW.x + static_cast<double>(i) * cellWidth;
		y[i] = cornerSW.y + static_cast<double>(i) * cellHeight;
	}
}

void Grid4x4::SetUniformZValue(double const zValue)
{
	z[0][0] = z[0][1] = z[0][2] = z[0][3] =
	z[1][0] = z[1][1] = z[1][2] = z[1][3] =
	z[2][0] = z[2][1] = z[2][2] = z[2][3] =
	z[3][0] = z[3][1] = z[3][2] = z[3][3] = zValue;
}

void Grid4x4::SetZValues(double const zValues[16])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			z[i][j] = *zValues;
			zValues++;
		}
	}
}

void Grid4x4::SetZValues(double const zValues[4][4])
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			z[i][j] = zValues[i][j];		
}

double Grid4x4::CellWidth() const
{
	return x[1] - x[0];
}

double Grid4x4::CellHeight() const
{
	return y[1] - y[0];
}

Vector2D Grid4x4::GridNode(unsigned int i, unsigned int j) const
{
	if (i > 3 || j > 3)
		return Vector2D();

	return Vector2D(x[i], y[j]);
}

Rect3D Grid4x4::GetCell(unsigned int i, unsigned int j) const
{
	if (i > 2 || j > 2)
		return Rect3D();

	Rect3D result(Rect(GridNode(i, j), CellWidth(), CellHeight()));
	result.z_SW = z[i][j];
	result.z_SE = z[i][j + 1];
	result.z_NW = z[i + 1][j];
	result.z_NE = z[i + 1][j + 1];

	return result;
}

#pragma endregion

#pragma region TimeSeries Definitions

TimeSeries::TimeSeries()
{
	size = 2;
	series = new std::pair<size_t, double>[size]();
}

TimeSeries::TimeSeries(size_t _size)
{
	size = Max(_size, static_cast<size_t>(2));
	series = new std::pair<size_t, double>[size](); 
}

TimeSeries::TimeSeries(std::vector<std::pair<size_t, double>> const & ts)
{
	if (ts.size() < 2)
	{
		size = 2;
		series = new std::pair<size_t, double>[size]();
		return;
	}

	size = ts.size();
	series = new std::pair<size_t, double>[size]();

	for (int i = 0; i < size; i++)
		series[i] = ts[i];
}

TimeSeries::~TimeSeries()
{
	if (series != NULL)
		delete[] series;
}

void TimeSeries::operator=(TimeSeries const & ts)
{
	delete[] series;
	
	size = ts.size;
	series = new std::pair<size_t, double>[size]();

	for (int i = 0; i < size; i++)
		series[i] = ts.series[i];
}

bool TimeSeries::IsValid() const
{
	//first element must be {0, 0.0}
	//no negative precipitation
	//every size_t is greater than preceding one.
	
	if (series[0].first != 0 || series[0].second != 0.0)
		return false;

	for (size_t i = 1; i < size; i++)
	{
		if (series[i].first <= series[i - 1].first
			|| series[i].second < 0.0)
			return false;
	}

	return true;
}

void TimeSeries::AdjustSize(size_t newSize)
{
	if (newSize == size) //nothing to do
		return;
	
	newSize = Max(newSize, static_cast<size_t>(2));

	//create a temporary holder for current data
	std::pair<size_t, double> * tempHolder = new std::pair<size_t, double>[newSize]();

	size_t minSize = Min(newSize, size);
	
	for (size_t i = 0; i < minSize; i++)
		tempHolder[i] = series[i];
	
	delete[] series;
	
	series = tempHolder;
	size = newSize;
}

double TimeSeries::HoursToLocalUnits(double time) const
{
	switch (timeUnit)
	{
	case TimeUnit::second:
		return time * 3600.0;
	case TimeUnit::minute:
		return time * 60.0;
	case TimeUnit::hour:
		return time;
	case TimeUnit::day:
		return time / 24.0;
	default:
		return 0.0;
	}
}

double TimeSeries::SampleRate(double timeSinceStart) const//, double timeSpan, InterpolationType interpolationType) const
{
	//Simplest way to convert incremental time-series to rate is to assume it constant between intervals, then divide by interval duration.
	//e.g., for a ts of (0, 0.0), (5, 10.0), (10, 15.0), (15, 5.0), the rate in first interval is (10 / (5 - 0)) = 2 mm/hr.
	//TODO Research better ways to address this.

	//convert timeSinceStart (in hours) to TS units
	double adjustedTimeSinceStart = HoursToLocalUnits(timeSinceStart);

	if (adjustedTimeSinceStart > series[size - 1].first || timeSinceStart < 0.0)
		return 0.0;

	size_t upperBound = size - 1;
	for (size_t i = 1; i < size; i++)
	{
		if (series[i].first >= adjustedTimeSinceStart)
		{
			upperBound = i;
			break;
		}
	}

	return HoursToLocalUnits(1.0) * series[upperBound].second / static_cast<double>(series[upperBound].first - series[upperBound - 1].first);

	//double relativePosition = (timeSinceStart - series[upperBound - 1].first) / (series[upperBound].first - series[upperBound - 1].first);
	//
	//switch (interpolationType)	
	//{
	//case InterpolationType::nearest:
	//	return relativePosition < 0.5 ? 0.0 : series[upperBound].second;
	//case InterpolationType::linear:
	//{
	//	double bounds[2];
	//	bounds[0] = 0.0;
	//	bounds[1] = series[upperBound].second;
	//	return LinearInterpolationNormalized(relativePosition, bounds);
	//}
	//case InterpolationType::cubic:
	//{
	//	double bounds[4];
	//	bounds[0] = bounds[1] = 0.0;
	//	bounds[2] = series[upperBound].second;
	//	bounds[3] = upperBound == size - 1 ? series[upperBound].second : series[upperBound + 1].second;
	//	return CubicInterpolationNormalized(relativePosition, bounds);
	//}
	//default: //shouldn't happen
	//	return 0.0;
	//}
}

#pragma endregion

#pragma region ModelParameters Definitions

ModelParameters::ModelParameters()
{
}

ModelParameters::~ModelParameters()
{
	/*if (unitTimeSeries != NULL)
		delete[] unitTimeSeries;*/
}

#pragma endregion


std::unique_ptr<double> ToUTM(double lng, double lat)
{
	//converting this http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html
	// to c++
	//also https://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.HTM
	//alternatively https://arxiv.org/pdf/1002.1417.pdf

	//std::cout << "\nDevWarning: Using Karney's method to convert from decimal degrees to UTM\n";

	unsigned int zone = (unsigned int)floor((lng + 180.0) / 6.0) + 1;
	double central_meridian_longitude = ((double)(zone - 1) * 6.0 - 180.0 + 3.0) * PI_CONSTANT / 180.0; //in radians
	//TODO consider Norway/Svalbard exceptions.

	const double e = sqrt(1.0 - pow((WGS_EARTH_RADIUS_POLES / WGS84_EARTH_RADIUS_EQUATOR), 2.0));
	const double n = (WGS84_ELIPSOID_FLATTENING / (2 - WGS84_ELIPSOID_FLATTENING));

	lat = lat * PI_CONSTANT / 180.0;
	lng = (lng * PI_CONSTANT / 180.0) - central_meridian_longitude;

	const double tao = tan(lat);
	const double sigma = sinh(e * atanh(e * tao / sqrt(1.0 + pow(tao, 2.0))));
	const double tao_prime = (tao * sqrt(1.0 + pow(sigma, 2.0))) - (sigma * sqrt(1.0 + pow(tao, 2.0)));
	const double xi_prime = atan(tao_prime / cos(lng));
	const double eta_prime = asinh(sin(lng) / sqrt(pow(tao_prime, 2.0) + pow(cos(lng), 2)));
	const double A = (WGS84_EARTH_RADIUS_EQUATOR / (1.0 + n)) * (1.0 + (1.0 / 4.0)*pow(n, 2.0) + (1.0 / 64.0)*pow(n, 4.0) + (1.0 / 256.0)*pow(n, 6.0)); //checked up to n^4 in source paper

	double alpha[6] =
	{ 
		(1.0 / 2.0) * n	- (2.0 / 3.0) * pow(n,2.0)		+ (5.0 / 16.0) * pow(n,3.0)		+ (41.0 / 180.0)*pow(n,4.0)			- (127.0 / 288.0)*pow(n,5.0)			+ (7891.0 / 37800.0)*pow(n, 6.0),
							  (13.0 / 48.0) * pow(n, 2.0)	- (3.0 / 5.0) * pow(n, 3.0)		+ (557.0 / 1440.0)*pow(n, 4.0)		+ (281.0 / 630.0)*pow(n, 5.0)		- (1983433.0 / 1935360.0)*pow(n, 6.0),
																  (61.0 / 240.0) * pow(n, 3.0)	- (103.0 / 140.0)*pow(n, 4.0)			+ (15061.0 / 26880.0)*pow(n, 5.0)	+ (167603.0 / 181440.0)*pow(n, 6.0),
																									  (49561.0 / 161280.0)*pow(n, 4.0)	- (179.0 / 168.0)*pow(n, 5.0)			+ (6601661.0 / 7257600.0)*pow(n, 6.0),
																																			  (34729.0 / 80640.0)*pow(n, 5.0)	- (3418889.0 / 1995840.0)*pow(n, 6.0),
																																													  (212378941.0 / 319334400.0)*pow(n, 6.0)
	};

	double xi = xi_prime;
	double eta = eta_prime;
	double p_prime = 1.0, q_prime = 0.0;
	for (int i = 0; i < 6; i++)
	{
		xi = xi + alpha[i] * sin(2.0 * (i + 1.0)*xi_prime) * cosh(2.0 * (i + 1.0) *eta_prime);
		eta = eta + alpha[i] * cos(2.0 * (i + 1.0) *xi_prime) * sinh(2.0 * (i + 1.0) *eta_prime);
		p_prime = p_prime + 2.0*i*cos(2.0 * (i + 1.0) *xi_prime) * cosh(2.0 * (i + 1.0)* eta_prime);
		q_prime = q_prime + 2.0*i*sin(2.0 * (i + 1.0)*xi_prime) * sinh(2.0 * (i + 1.0) * eta_prime);
	}

	double x = UTM_MERIDIAN_SCALE * A * eta;
	double y = UTM_MERIDIAN_SCALE * A * xi;

	x = x + UTM_FALSE_EASTING;
	if (y < 0) y = y + UTM_FALSE_NORTHING; //in case the point was in sourthern hemisphere. for norther hemi, the y above is ok.

	std::unique_ptr<double> coords = std::unique_ptr<double>(new double[2]);
	coords.get()[0] = y;
	coords.get()[1] = x;

	//The part bellow would severely impact performance, in case of large profiles.
	/*if (isDebug)
		std::cout << "\n in ToUTM, returning coords: " << coords.get()[0] << " and " << coords.get()[1];*/

	return coords;
}

std::unique_ptr<double> ToWGS84(double easting, double northing, bool isNortherHemisphere, int zone)
{
	//Using the some source-code referenced in ToUTM() and converting it to C++

	double _easting = easting - UTM_FALSE_EASTING;
	double _northing = isNortherHemisphere ? northing : northing - UTM_FALSE_NORTHING;

	double eccentricity = sqrt(WGS84_ELIPSOID_FLATTENING * (2.0 - WGS84_ELIPSOID_FLATTENING)); //this could be calculted outside and hardcoded into this program.
	double n = WGS84_ELIPSOID_FLATTENING / (2.0 - WGS84_ELIPSOID_FLATTENING); //ditto

	double n2 = n * n;
	double n3 = n2 * n;
	double n4 = n3 * n;
	double n5 = n4 * n;
	double n6 = n5 * n;

	double A = (WGS84_EARTH_RADIUS_EQUATOR / (1.0 + n)) * (1.0 + n2 * (1.0 / 4.0) + n4 * (1.0 / 64.0) + n6 * (1.0 / 256.0));

	double eta = _easting / (UTM_MERIDIAN_SCALE * A);
	double xi = _northing / (UTM_MERIDIAN_SCALE * A);

	double beta[6] =
	{
		 n * (1.0 / 2.0)	- n2 * (2.0 / 3.0)	+ n3 * (37.0 / 96.0)	- n4 * (1.0 / 360.0)		- n5 * (81.0 / 512.0)		+ n6 * (96199.0 / 604800.0),
							  n2 * (1.0 / 48.0)	+ n3 * (1.0 / 15.0)		- n4 * (437.0 / 1440.0)		+ n5 * (46.0 / 105.0)		- n6 * (1118711.0 / 3870720.0),
												  n3 * (17.0 / 480.0)	- n4 * (37.0 / 840.0)		- n5 * (209.0 / 4480.0)		+ n6 * (5569.0 / 90720.0),
																		  n4 * (4397.0 / 161280.0)	- n5 * (11.0 / 504.0)		- n6 * (830251.0 / 7257600.0),
																									  n5 * (4583.0 / 161280.0)	- n6 * (108847.0 / 3991680.0),
																																  n6 * (20648693.0 / 638668800.0)
	};

	double xi_prime = xi;
	double eta_prime = eta;

	for (int i = 0; i < 6; i++)
	{
		xi_prime -= beta[i] * sin(2.0 * (i + 1) * xi) * cosh(2.0 * (i + 1) * eta);
		eta_prime -= beta[i] * cos(2.0 * (i + 1) * xi) * sinh(2.0 * (i + 1) * eta);
	}

	double sinh_eta_prime = sinh(eta_prime);
	double sin_xi_prime = sin(xi_prime);
	double cos_xi_prime = cos(xi_prime);

	double tau_prime = sin_xi_prime / sqrt(sinh_eta_prime * sinh_eta_prime + cos_xi_prime * cos_xi_prime);

	double deltaTau = 0.1;
	double tau = tau_prime;

	while (abs(deltaTau) > 0.00000000001)
	{
		double sigma = sinh(eccentricity * atanh(eccentricity * tau / sqrt(1.0 + tau * tau)));
		double tau_i_prime = tau * sqrt(1.0 + sigma * sigma) - sigma * sqrt(1 + tau * tau);
		deltaTau = ((tau_prime - tau_i_prime) / sqrt(1.0 + tau_i_prime * tau_i_prime)) * ((1.0 + (1.0 - eccentricity * eccentricity) * tau * tau) / ((1.0 - eccentricity * eccentricity) * sqrt(1.0 + tau * tau)));
		tau += deltaTau;
	}

	double phi = atan(tau);
	double lambda = atan2(sinh_eta_prime, cos_xi_prime);

	double p = 1.0;
	double q = 0.0;

	for (int i = 0; i < 6; i++)
	{
		p -= 2.0 * i * beta[i] * cos(2.0 * (i + 1) * xi) * cosh(2.0 * (i + 1) * eta);
		q += 2.0 * i * beta[i] * sin(2.0 * (i + 1) * xi) * sinh(2.0 * (i + 1) * eta);
	}

	double gamma_prime = atan(tan(xi_prime) * tanh(eta_prime));
	double gamme_prime_prime = atan2(q, p);
	double gamma = gamma_prime + gamme_prime_prime;

	double sin_phi = sin(phi);

	double k_prime = sqrt(1.0 - eccentricity * eccentricity *sin_phi * sin_phi) * sqrt(1.0 + tau * tau) * sqrt(sinh_eta_prime * sinh_eta_prime + cos_xi_prime * cos_xi_prime);
	double k_prime_prime = (A / WGS84_EARTH_RADIUS_EQUATOR) * sqrt(p * p + q * q);
	double k = UTM_MERIDIAN_SCALE * k_prime * k_prime_prime;

	double lambda_0 = (PI_CONSTANT / 180.0)* (double)((zone - 1) * 6 - 180 + 3);
	lambda += lambda_0;

	std::unique_ptr<double> coords = std::unique_ptr<double>(new double[2]);

	coords.get()[0] = lambda * (180.0 / PI_CONSTANT);
	coords.get()[1] = phi * (180.0 / PI_CONSTANT);


	double convergence = gamma * 180.0 / PI_CONSTANT;
	double scale = k;

	return coords;
}

Vector2D ProjectPoint(Vector2D const & point)
{
	std::unique_ptr<double> result = ToUTM(point.x, point.y);
	return Vector2D(result.get()[0], result.get()[1]);
}

Vector2D WrapPoint(Vector2D const & point, bool isNorthernHemisphere, int zone)
{
	std::unique_ptr<double> result = ToWGS84(point.x, point.y, isNorthernHemisphere, zone);
	return Vector2D(result.get()[0], result.get()[1]);
}

double LinearInterpolationNormalized(double normalizedPoint, double const values[2])
{
	return normalizedPoint * values[1] + (1.0 - normalizedPoint) * values[0];
}

double CubicInterpolationNormalized(double const normalizedPoint, double const values[4])
{
	return values[1] + 0.5 * normalizedPoint *
			(values[2] - values[0] + normalizedPoint *
				(2.0 * values[0] - 5.0 * values[1] + 4.0 * values[2] - values[3] + normalizedPoint *
					(3.0 * (values[1] - values[2]) + values[3] - values[0])));
}

double BilinearInterpolation(Vector2D const & point, Grid4x4 const & grid)
{
	return BilinearInterpolation(point, grid.GetCell(1, 1));
}

double BilinearInterpolation(Vector2D const & point, Rect3D const & grid)
{
	Vector2D normalizedPoint = point.Normalize(grid.MinCorner(), grid.MaxCorner());
	double linearIntX[2];
	
	double edgeS[2]{ grid.z_SW, grid.z_SE };
	double edgeN[2]{ grid.z_NW, grid.z_NE };

	linearIntX[0] = LinearInterpolationNormalized(normalizedPoint.x, edgeS);
	linearIntX[1] = LinearInterpolationNormalized(normalizedPoint.x, edgeN);

	return LinearInterpolationNormalized(normalizedPoint.y, linearIntX);
}

double BicubicInterpolation(Vector2D const & point, Grid4x4 const & grid)
{
	//reference:
	//http://www.paulinternet.nl/?page=bicubic

	Vector2D normalizedPoint = point.Normalize(grid.GridNode(1, 1), grid.GridNode(2, 2));
	double cubicIntX[4];

	for (int i = 0; i < 4; i++)
		cubicIntX[i] = CubicInterpolationNormalized(normalizedPoint.x, grid.z[i]);

	return CubicInterpolationNormalized(normalizedPoint.y, cubicIntX);
}

