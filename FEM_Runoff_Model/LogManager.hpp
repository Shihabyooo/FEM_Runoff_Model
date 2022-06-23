#pragma once
//#include <wchar.h>
//#include <vector>
//#include <filesystem>
//#include <iostream>
#include <string>

#include "Globals.hpp"

//#include "GUI_Requirements.hpp"

void DrawLogPane();
namespace LogMan
{
	void Init(bool logToDisk = true);
	void Terminate();
	void Log(std::string const & content, LogEntryType type = LOG_NORM);
	void Log(char const * content, LogEntryType type = LOG_NORM);
}