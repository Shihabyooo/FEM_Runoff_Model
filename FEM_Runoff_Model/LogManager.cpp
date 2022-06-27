#include <chrono>
#include <ctime>

#include "LogManager.hpp"
#include "GUI_Requirements.hpp"
#include "FileIO.hpp"

std::vector<LogEntry> logHistory;
bool canLogToDisk = false;

void DrawLogEntry(LogEntry const * entry)
{
	switch (entry->type)
	{
	case LogEntryType::normal:
		//ImGui::TextColored(ImVec4(0.0f, 0.0f, 0.0f, 1.0f), logHistory[order].content.c_str());
		//ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
		ImGui::TextWrapped(entry->content.c_str());
		//ImGui::PopStyleColor();
		return;
	case LogEntryType::warning:
		//ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), logHistory[order].content.c_str());
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.5f, 0.0f, 1.0f));
		ImGui::TextWrapped(entry->content.c_str());
		ImGui::PopStyleColor();
		return;
	case LogEntryType::error:
		//ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), logHistory[order].content.c_str());
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
		ImGui::TextWrapped(entry->content.c_str());
		ImGui::PopStyleColor();
		return;
	case LogEntryType::success:
		//ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), logHistory[order].content.c_str());
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 1.0f, 0.0f, 1.0f));
		ImGui::TextWrapped(entry->content.c_str());
		ImGui::PopStyleColor();
		return;
	default:
		//ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), logHistory[order].content.c_str());
		ImGui::TextWrapped(entry->content.c_str());
		return;
	}

}

void DrawLogPane()
{
	ImGuiWindowFlags windowFlags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	ImGui::SetNextWindowPos(ImVec2(logPaneDimensions.positionX, logPaneDimensions.positionY), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(logPaneDimensions.width, logPaneDimensions.height), ImGuiCond_Always);

	ImGui::Begin("Log", NULL, windowFlags);

	ImGui::PushItemWidth(ImGui::GetFontSize() * -12); //Use fixed width for labels (by passing a negative value), the rest goes to widgets. We choose a width proportional to our font s

	//ImGui::TextWrapped("Log goes here");

	for (auto it = logHistory.begin(); it != logHistory.end(); ++it)
		DrawLogEntry(&(*it));

	ImGui::End();
}

std::string GetTimeStamp()
{
	std::time_t timeNow;
	timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

	tm tmStruct;// = gmtime(&timeNow);
	gmtime_s(&tmStruct, &timeNow);

	std::string timeSamp = std::to_string(tmStruct.tm_year + 1900) + "/" + std::to_string(tmStruct.tm_mon + 1) + "/" + std::to_string(tmStruct.tm_mday) + " "
		+ std::to_string(tmStruct.tm_hour) + ":" + std::to_string(tmStruct.tm_min) + ":" + std::to_string(tmStruct.tm_sec) + " | ";
	return timeSamp;
}

void LogMan::Init(bool logToDisk)
{
	canLogToDisk = logToDisk && FileIO::InitLogFile();
}

void LogMan::Terminate()
{
	FileIO::CloseLogFile();
}

void LogMan::Log(std::string const & content, LogEntryType type)
{
	//TODO enrich the logging mechanism with a max value, after which the oldest logs are removed.
	logHistory.push_back(LogEntry(GetTimeStamp() + content, type));
	
	if(canLogToDisk)
		FileIO::WriteToLog(logHistory.back());

#ifdef LOG_TO_CLI
	std::cout << content.c_str() << std::endl;
#endif
}

void LogMan::Log(char const * content, LogEntryType type)
{
	Log(std::string(content), type);
}
