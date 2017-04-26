#include <iostream>
#include <fstream>

#include "json.hpp"

int main(int argc, char *argv[])
{
	// Read in arguments
	if (argc < 2)
	{
		std::cout << "Need input param file. Exiting now..." << std::endl;
		return 0;
	}

	// Read json file in memory
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);
	jsonFile >> jsonInput;


	// Output
	std::cout << jsonInput["program-info"]["debug-log-path"] << std::endl;


	return 0;
}

// Run:
// ./ScalingFramework ../input/input.json