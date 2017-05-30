#ifndef _RANDOM_UTILS_H_
#define _RANDOM_UTILS_H_


#include <string>

#ifdef __linux__ 
	#include <sys/resource.h>
	#include <sys/time.h>
#elif _WIN32

#elif __APPLE__

#endif



std::string toBinary(int num)
{
    std::string binStr = "";

    while (num>0)
	{
	    binStr = std::to_string(num % 2) + binStr;
	    num /= 2;
	}

	return binStr;
}


int toDecimal(std::string binary)
{
	int num = 0;
	int pow = 1;

	for (int i=binary.size()-1; i>=0; i--)
	{
		num += ( (binary[i] == '0') ? 0 : 1 ) * pow;
		pow *= 2;
	}

	return num;
}


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

