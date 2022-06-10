#pragma once
#define PROGRAM_NAME "FEM_Runoff_Model"

//error codes
#define SUCCESS 0
#define FAILED_MAIN_WINDOW_INITIALIZE 1
#define FAILED_VIEWPORT_CREATE 2

//struct Vector3
//{
//public:
//	Vector3(float _x, float _y, float _z)
//	{
//		x = _x;
//		y = _y;
//		z = _z;
//	}
//
//	float x;
//	float y;
//	float z;
//};

struct Vector2
{
public:
	Vector2()
	{
		x = 0;
		y = 0;
	};

	Vector2(float _x, float _y)
	{
		x = _x;
		y = _y;
	};
	
	Vector2 operator+(Vector2 const & vec2)
	{
		return (Vector2(x + vec2.x, y + vec2.y));
	}

	Vector2 operator-(Vector2 const & vec2)
	{
		return (Vector2(x - vec2.x, y - vec2.y));
	};

	Vector2 operator*(float const & scalar)
	{
		return (Vector2(x * scalar, y * scalar));
	};

	float x;
	float y;
};

static float Min(float  const & a, float  const & b)
{
	//std::cout << "comparing least of " << a << " and " << b << " result: " << (a > b ? b : a) << std::endl;
	return (a > b ? b : a);
}

static float Max(float  const & a, float  const & b)
{
	//std::cout << "comparing greater of " << a << " and " << b << " result: " << (a > b ? a : b) << std::endl;
	return (a > b ? a : b);
}

static float Clamp(float const & a, float  const & b, float  const & c)
{
	float min = Min(b, c);
	float max = Max(b, c);
	return (a > max ? max : (a < min ? min : a));
}

static void Print(Vector2 & vec)
{
	std::cout << vec.x << ", " << vec.y << "\n";
}