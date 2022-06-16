#include "Main.hpp"
#include "LogManager.hpp"

int main(int argc, char ** argv)
{
	LogMan::Log("Startup.");
	return StartUI(1280, 720);
}