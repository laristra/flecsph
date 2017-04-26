#ifndef _RANDOM_UTILS_H_
#define _RANDOM_UTILS_H_


#include <string>

#ifdef __linux__ 
	#include <sys/resource.h>
	#include <sys/time.h>
#elif _WIN32

#elif __APPLE__

#endif




#ifdef __linux__ 

inline size_t getCurrentMemoryUse()
{
	// http://www.tutorialspoint.com/unix_system_calls/getrusage.htm
	// http://appcrawler.com/wordpress/2013/05/13/simple-example-of-tracking-memory-using-getrusage/
	struct rusage _rusage;
	getrusage(RUSAGE_SELF,&_rusage);

	return (size_t) (_rusage.ru_maxrss);
}

#elif _WIN32

#else

#endif





#endif

