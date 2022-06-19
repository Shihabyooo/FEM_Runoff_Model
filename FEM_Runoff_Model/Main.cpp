#include "Main.hpp"
#include "LogManager.hpp"

int main(int argc, char ** argv)
{
	LogMan::Log("Startup.");
	return StartUI(1280, 720);

	////testing inteprolation
	//double width = 2.0f;
	//double height = 2.0f;
	//Grid4x4 grid(Vector2D(5.0f, 5.0f), width, height); //from 5.0 to 9.0 on each axes
	//double testZ[4][4] = {	{1.0f, 1.0f, 1.0f, 1.0f},
	//						{1.0f, 2.0f, 2.0f, 1.0f},
	//						{1.0f, 2.0f, 2.0f, 1.0f},
	//						{1.0f, 1.0f, 1.0f, 1.0f} };
	//
	//grid.SetZValues(testZ);



	//double x, y;
	//int subdivs = 5;
	//for (int i = 0; i <= subdivs; i++)
	//{
	//	x = 5.0f + width + ((double)(subdivs - i) / (double)subdivs) * width;
	//	y = 5.0f + height + ((double)(subdivs - i) / (double)subdivs) * height;

	//	std::cout << "blinear: " << x << ", " << y << ": " << BilinearInterpolation(Vector2D(x, y), grid) << std::endl;
	//	std::cout << "bicubic: " << x << ", " << y << ": " << BicubicInterpolation(Vector2D(x, y), grid) << std::endl;
	//	std::cout << "\n";
	//}
	//std::cin >> x;
}