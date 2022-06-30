#include "Main.hpp"
#include "LogManager.hpp"

#include "Solvers.hpp" //only needed for testing
#include <chrono> //ditto

#define NOW std::chrono::high_resolution_clock::now()
#define TIME std::chrono::high_resolution_clock::time_point
#define DURATION_SINCE(x) std::chrono::duration_cast<std::chrono::milliseconds>(NOW - x).count()

//Ax = b
Matrix_f64 testA;
Vector_f64 testB;
Vector_f64 testR, testX;
int size = 100;
int modulus = 100;

void GenerateTestSystem(bool isSymmetric, bool isDiagDominant)
{
	testA = Matrix_f64(size, size);
	testB = Vector_f64(size);
	//testR and testX are set by solver.
	srand(time(0));
	double tempVal1;
	
	tempVal1 = testA[0][1] = rand() % modulus;
	testA[0][0] = testA[0][1] + 1;
	/*testA[0][0] = 2.5f;
	testA[0][1] = 1.0f;*/

	testB[0] = rand() % modulus;
	
	for (int i = 1; i < size - 1; i++)
	{
		testB[i] = rand() % modulus;

		if (isSymmetric)
			testA[i][i - 1] = tempVal1;
		else
			//testA[i][i - 1] = rand() % modulus;
			testA[i][i - 1] = 1.6f;
		
		//tempVal1 = testA[i][i + 1] = rand() % modulus;
		tempVal1 = testA[i][i + 1] = 10.0f;

		if (isDiagDominant)
			testA[i][i] = testA[i][i - 1] + testA[i][i + 1] + rand() % modulus;
		else
			testA[i][i] = rand() % modulus;

	
		/*testA[i][i - 1] = 1.0f;
		testA[i][i] = 2.5f;
		testA[i][i + 1] = 1.0f;*/
	}

	testA[size - 1][size - 2] = tempVal1;
	testA[size - 1][size - 1] = rand() % modulus;
	/*testA[size - 1][size - 2] = 1.0f;
	testA[size - 1][size - 1] = 2.5f;*/
	
	testB[size - 1] = rand() % modulus;

	//testA.DisplayArrayInCLI(1);
	//testB.DisplayArrayInCLI(1);
}

void SolverTest()
{
	TIME start;
	long time;

	for (int i = 0; i < 10; i++)
	{
		std::cout << "\n================================================\n";
		std::cout << "Test iteration: " << i << std::endl;
		std::cout << "\================================================\n";
		GenerateTestSystem(true, true);

		//start = NOW;
		//SolverSimple(testA, testB, testX, testR);
		/*time = DURATION_SINCE(start);
		std::cout << "Residual: " << testR.Magnitude() << "\tDuration: " << time << "ms" << std::endl << std::endl;*/

		/*start = NOW;
		SolverJacobi(testA, testB, testX, testR);
		time = DURATION_SINCE(start);
		std::cout << "Residual: " << testR.Magnitude() << "\tDuration: " << time << "ms" << std::endl << std::endl;

		start = NOW;
		SolverSOR(testA, testB, testX, testR);
		time = DURATION_SINCE(start);
		std::cout << "Residual: " << testR.Magnitude() << "\tDuration: " << time << "ms" << std::endl << std::endl;*/

		start = NOW;
		SolverPCG(testA, testB, testX, testR);
		time = DURATION_SINCE(start);
		std::cout << "Residual: " << testR.Magnitude() << "\tDuration: " << time << "ms" << std::endl << std::endl;

		start = NOW;
		SolverBiCG(testA, testB, testX, testR);
		time = DURATION_SINCE(start);
		std::cout << "Residual: " << testR.Magnitude() << "\tDuration: " << time << "ms" << std::endl << std::endl;

		start = NOW;
		SolverCGS(testA, testB, testX, testR);
		time = DURATION_SINCE(start);
		std::cout << "Residual: " << testR.Magnitude() << "\tDuration: " << time << "ms" << std::endl << std::endl;
	}
}

int main(int argc, char ** argv)
{
	LogMan::Init(true);
	LogMan::Log("Startup.");
	int returnVal = StartUI(1280, 720);
	
	//SolverTest();

	LogMan::Terminate();
	//std::cin.sync(); //TODO remove
	//std::cin.get(); //TODO remove

	return returnVal;
}
