#pragma once
#include <iostream>
#define PROGRAM_NAME "FEM_Runoff_Model"

//error codes
#define SUCCESS 0
#define FAILED_MAIN_WINDOW_INITIALIZE 1
#define FAILED_VIEWPORT_CREATE 2

//Log Types definitions
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

struct Vector2
{
public:
	Vector2()
	{
		x = 0.0f;
		y = 0.0f;
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

	Vector2 & operator= (Vector2 const & vec2)
	{
		x = vec2.x;
		y = vec2.y;
		return *this;
	}

	float x;
	float y;
};

struct Vector2Int
{
public:
	Vector2Int()
	{
		x = 0;
		y = 0;
	};

	Vector2Int(int _x, int _y)
	{
		x = _x;
		y = _y;
	};

	Vector2Int operator+(Vector2 const & vec2)
	{
		return (Vector2Int(x + vec2.x, y + vec2.y));
	}

	Vector2Int operator-(Vector2 const & vec2)
	{
		return (Vector2Int(x - vec2.x, y - vec2.y));
	};

	Vector2Int operator*(int const & scalar)
	{
		return (Vector2Int(x * scalar, y * scalar));
	};

	int x;
	int y;
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
static float Min(float const & a, float const & b)
{
	return (a > b ? b : a);
}

static int Min(int const & a, int const & b)
{
	return (a > b ? b : a);
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