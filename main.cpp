#include <iostream>
#include <fstream>

#include "json.hpp"
#include "octree.h"
#include "scalingTest.h"


int main(int argc, char *argv[])
{
	// Read in arguments
	if (argc < 2)
	{
		std::cout << "Need input json file. Exiting now..." << std::endl;
		return 0;
	}


	// Read json file in memory
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);
	jsonFile >> jsonInput;


	//// Output
	//std::cout << jsonInput["program-info"]["debug-log-path"] << std::endl;

	Octree testOctree;
	ScalingTest ioTesting;

	testOctree.buildTree( jsonInput["data"]["num-octree-levels"] );
	std::cout << testOctree.print();

	MPI_Init(NULL, NULL);
	
	ioTesting.initMPIScaling(MPI_COMM_WORLD);
	ioTesting.runScalingTest(jsonInput["data"]["num-particles"], jsonInput["data"]["num-timesteps"], testOctree, jsonInput["output"]["filename"]);
	std::cout << ioTesting.getTimingLog() << std::endl;
	

	MPI_Finalize();

	return 0;
}

// Run:
// ./ScalingFramework ../input/input.json