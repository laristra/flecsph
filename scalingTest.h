#ifndef _SCALING_TESTS_H_
#define _SCALING_TESTS_H_


#include <string>
#include <iostream>
#include <random>

#include <mpi.h>

#include "timer.h"
#include "log.h"
#include "hdf5SimIO.h"

class ScalingTest
{
	int myRank;
	int numRanks;
	bool theadingON;
	MPI_Comm mpiComm;
	Log timingLog;

  public:
  	ScalingTest(){ theadingON=false; myRank=0; numRanks=1; }
  	~ScalingTest(){};

  	void enableTheading(){ theadingON=true; }
  	void initMPIScaling(MPI_Comm _comm);
  	void runScalingTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename);

  	std::string getTimingLog(){ return timingLog.getLog(); }
};


inline void ScalingTest::initMPIScaling(MPI_Comm _comm)
{
	mpiComm = _comm;
	MPI_Comm_rank(_comm, &myRank);
	MPI_Comm_size(_comm, &numRanks);

	std::cout << myRank << std::endl;
}


inline void ScalingTest::runScalingTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename)
{
	Timer dataCreation, dataWriting;

	int myNumParticles = _numParticles/numRanks;

	//
	// Get Partition from octree
	std::string nodeID = simOct.getNodeAt( (size_t)myRank );

	float nodeExtents[6];
	simOct.getSpatialExtent(nodeID, nodeExtents);
	std::cout << myRank << " ~ # particles: " << myNumParticles << ", Extents: "	
										<< nodeExtents[0] << " - " << nodeExtents[1] << "   "
										<< nodeExtents[2] << " - " << nodeExtents[2] << "   "
										<< nodeExtents[4] << " - " << nodeExtents[5] << std::endl;


	for (int t=0; t<numTimesteps; t++)
	{
		//
		// Generate Data
	  dataCreation.start();
		// Some variables
		float *_x_data = new float[myNumParticles];
		float *_y_data = new float[myNumParticles];
		float *_z_data = new float[myNumParticles];
		float *_pressure_data = new float[myNumParticles];


		std::default_random_engine eng( (std::random_device()) () );
		std::uniform_real_distribution<float> xDis(nodeExtents[0], nodeExtents[1]);
		std::uniform_real_distribution<float> yDis(nodeExtents[2], nodeExtents[3]);
		std::uniform_real_distribution<float> zDis(nodeExtents[4], nodeExtents[5]);
		std::uniform_real_distribution<float> pressureDis(-1000.0, 1000.0);
		for (size_t i=0; i<myNumParticles; i++)
		{
			_x_data[i] = xDis(eng);
			_y_data[i] = yDis(eng);
			_z_data[i] = zDis(eng);
			_pressure_data[i] = pressureDis(eng);
		}
	  dataCreation.stop();


	  dataWriting.start();
		//
		// Create hdf5
		Flecsi_Sim_IO::HDF5SimIO testDataSet( filename.c_str() );

		testDataSet.createDataset(Flecsi_Sim_IO::unStructuredGrid, 4, 3);	// 4: num vars, 3: num dims

		Flecsi_Sim_IO::Variable _x, _y, _z, _pressure;
		_x.createVariable("x", Flecsi_Sim_IO::point, "float", myNumParticles, _x_data);	
		_y.createVariable("y", Flecsi_Sim_IO::point, "float", myNumParticles, _y_data);
		_z.createVariable("z", Flecsi_Sim_IO::point, "float", myNumParticles, _z_data);
		_pressure.createVariable("pressure", Flecsi_Sim_IO::point, "float", myNumParticles, _pressure_data);

		testDataSet.vars.push_back(_x);
		testDataSet.vars.push_back(_y);
		testDataSet.vars.push_back(_z);
		testDataSet.vars.push_back(_pressure);
	  
		testDataSet.writePointData(t, mpiComm);
	  dataWriting.stop();

	  	timingLog.addLog("Data creation for timestep " + std::to_string(t) + " took: " + std::to_string( dataCreation.getDuration() ) + " s.\n");
	  	timingLog.addLog("Data write for timestep " + std::to_string(t)+ " took: " + std::to_string( dataWriting.getDuration() ) + " s.\n\n");

		// Delete 
		if (_x_data != NULL)
			delete [] _x_data;
		if (_y_data != NULL)
			delete [] _y_data;
		if (_z_data != NULL)
			delete [] _z_data;
		if (_pressure_data != NULL)
			delete [] _pressure_data;

		// Sync writing before going to the next step
		MPI_Barrier(mpiComm);
	}

}

#endif