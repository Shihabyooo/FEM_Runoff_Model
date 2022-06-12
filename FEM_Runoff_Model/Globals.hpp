#pragma once
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

static void Print(Vector2 & vec)
{
	std::cout << vec.x << ", " << vec.y << "\n";
}