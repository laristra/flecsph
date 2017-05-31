#include <iostream>
#include <fstream>

#include "json.hpp"
#include "octree.h"
#include "scalingTest.h"
#include "log.h"


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


	Octree testOctree;
	ScalingTest ioTesting;

	testOctree.buildTree( jsonInput["data"]["num-octree-levels"] );


	MPI_Init(NULL, NULL);
	int myRank, numRanks;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

	Log debugLog( std::to_string(myRank) + "_of_" + std::to_string(numRanks) + ".log" );


	ioTesting.initMPIScaling( MPI_COMM_WORLD );
	ioTesting.setIterationCount( jsonInput["output"]["repeats"] );
	ioTesting.runScalingTest( jsonInput["data"]["num-particles"], jsonInput["data"]["num-timesteps"], testOctree, jsonInput["output"]["filename"] );

	debugLog.addLog( ioTesting.getTimingLog() );
	debugLog.writeLogToDisk();

	

	MPI_Finalize();

	return 0;
}

// Run:
// mpirun -np 8 ./ScalingFramework ../input/input.json
