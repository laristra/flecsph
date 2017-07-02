#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>
#include <string>
#include <ctime>
#include <sstream>

class Timer
{
	std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
	std::chrono::duration<double> elapsed_seconds;

  public:

	Timer();
	~Timer();

	void start();
	void stop();
	double getDuration();					// time in seconds

	static std::string getCurrentTime();	// get the current time
};

inline Timer::Timer(){}
inline Timer::~Timer(){}

inline void Timer::start(){ startTime = std::chrono::system_clock::now(); }
inline void Timer::stop(){ endTime = std::chrono::system_clock::now(); elapsed_seconds = endTime - startTime; }
inline double Timer::getDuration(){ return elapsed_seconds.count(); }


inline std::string Timer::getCurrentTime()
{
	time_t now = time(0);
	tm *ltm = localtime(&now);

	std::stringstream ss;
	ss << "_"<< 1 + ltm->tm_mon << "_" << ltm->tm_mday << "__" << ltm->tm_hour << "_" << ltm->tm_min << "_" << ltm->tm_sec << "_" << std::endl;
	return ss.str();
}

#endif
