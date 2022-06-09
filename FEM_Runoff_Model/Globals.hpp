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
	
	Vector2 operator-(Vector2 vec2)
	{
		return (Vector2(x - vec2.x, y - vec2.y));
	};

	float x;
	float y;
};

static float Min(float a, float b)
{
	return a > b ? b : a;
}

static float Max(float a, float b)
{
	return a > b ? a : b;
}