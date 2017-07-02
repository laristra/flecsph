#ifndef _LOG_H_
#define _LOG_H_

#include <fstream>
#include <string>
#include <sstream>


class Log
{
	std::string outputFilename;
	std::string logMsg;

  public:
	Log();
	Log(std::string _outputFilename): outputFilename(_outputFilename){ logMsg = ""; }
	~Log();

	void setOutputFilename(std::string _outputFilename);
	void clearLog();
	void addLog(std::string _msg);
	void writeLogToDisk();
	std::string getLog();
};



inline Log::Log()
{
	outputFilename = "untitled.log";
	logMsg = "";
}

inline Log::~Log()
{
	outputFilename = "";
	logMsg = "";
}


inline void Log::setOutputFilename(std::string _outputFilename)
{
	outputFilename = _outputFilename;
}


inline void Log::clearLog()
{
	logMsg = "";
}

inline void Log::addLog(std::string _msg)
{
	logMsg += _msg;
}


inline void Log::writeLogToDisk()
{
  #ifdef ENABLE_DEBUG
	std::ofstream outputFile( outputFilename.c_str(), std::ios::out);
	outputFile << logMsg;
	outputFile.close();
  #endif //ENABLE_DEBUG
}


inline std::string Log::getLog()
{
	return logMsg;
}

#endif