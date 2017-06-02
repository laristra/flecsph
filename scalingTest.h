#ifndef _SCALING_TESTS_H_
#define _SCALING_TESTS_H_


#include <string>
#include <iostream>
#include <sstream>
#include <random>

#include <mpi.h>

#include "timer.h"
#include "log.h"
#include "hdf5ParticleIO.h"


class ScalingTest
{
	int myRank;
	int numRanks;
	bool theadingON;
	int iterationCount;
	MPI_Comm mpiComm;
	std::stringstream logStream;
	//Flecsi_Sim_IO::HDF5SimIO testDataSet;
	Flecsi_Sim_IO::HDF5ParticleIO testDataSet;

  public:
  	ScalingTest(){ theadingON=false; myRank=0; numRanks=1; iterationCount=1; }
  	~ScalingTest(){};

  	void enableTheading(){ theadingON=true; }
  	void setIterationCount(int _iterationCount){ iterationCount=_iterationCount; }
  	void initMPIScaling(MPI_Comm _comm);

  	void runScalingTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename);
  	void createPseudoData(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename);
  	void runWriteTest(int timeStep, size_t _numParticles);

  	std::string getTimingLog(){ return logStream.str(); }
};



inline void ScalingTest::initMPIScaling(MPI_Comm _comm)
{
	mpiComm = _comm;
	MPI_Comm_rank(_comm, &myRank);
	MPI_Comm_size(_comm, &numRanks);
}


inline void ScalingTest::runScalingTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename)
{
	createPseudoData(_numParticles, numTimesteps, simOct, filename);
}


inline void ScalingTest::createPseudoData(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename)
{
	Timer dataCreation, dataConversion;

	int myNumParticles = _numParticles/numRanks;

	//
	// Get Partition from octree
	std::string nodeID = simOct.getNodeAt( (size_t)myRank );


	float nodeExtents[6];
	simOct.getSpatialExtent(nodeID, nodeExtents);

	testDataSet.createDataset( filename.c_str(), mpiComm );

	testDataSet.writeGlobalAttribute("numParticles", "int64_t", 80);
	testDataSet.writeGlobalAttribute("gravitational constant", "double", 6.67300E-11);
	testDataSet.writeGlobalAttribute("numDims", "int32_t", 3);
	testDataSet.writeGlobalAttribute("use_fixed_timestep", "int32_t", 0);

	char *simName = "random_data";
	testDataSet.writeGlobalAttributeArray("sim_name", "string", simName);
	testDataSet.closeFile();



	for (int t=0; t<numTimesteps; t++)
	{
		//
		// Generate Data
	  dataCreation.start();
		// Some variables
		float *_x_data = new float[myNumParticles];
		float *_y_data = new float[myNumParticles];
		float *_z_data = new float[myNumParticles];
		double *_pressure_data = new double[myNumParticles];

		std::default_random_engine eng( (std::random_device()) () );
		std::uniform_real_distribution<float> xDis(nodeExtents[0], nodeExtents[1]);
		std::uniform_real_distribution<float> yDis(nodeExtents[2], nodeExtents[3]);
		std::uniform_real_distribution<float> zDis(nodeExtents[4], nodeExtents[5]);
		std::uniform_real_distribution<double> pressureDis(-1000.0, 1000.0);
		for (size_t i=0; i<myNumParticles; i++)
		{
			_x_data[i] = xDis(eng);
			_y_data[i] = yDis(eng);
			_z_data[i] = zDis(eng);
			_pressure_data[i] = pressureDis(eng);
		}
	  dataCreation.stop();



	  	//
		// Create hdf5
	  dataConversion.start();

	  	testDataSet.openFile(mpiComm);
	  	testDataSet.setTimeStep(t);


	  	Flecsi_Sim_IO::Attribute timeValue;
	  	timeValue.createAttribute("time", Flecsi_Sim_IO::timestep, "float", (float)t/10.0);
	  	testDataSet.timestepAttributes.push_back(timeValue);

		Flecsi_Sim_IO::Variable _x, _y, _z, _pressure;
		_x.createVariable("x", Flecsi_Sim_IO::point, "float", myNumParticles, _x_data);	
		_y.createVariable("y", Flecsi_Sim_IO::point, "float", myNumParticles, _y_data);
		_z.createVariable("z", Flecsi_Sim_IO::point, "float", myNumParticles, _z_data);
		_pressure.createVariable("pressure", Flecsi_Sim_IO::point, "double", myNumParticles, _pressure_data);

		testDataSet.vars.push_back(_x);
		testDataSet.vars.push_back(_y);
		testDataSet.vars.push_back(_z);
		testDataSet.vars.push_back(_pressure);
	  dataConversion.stop();

	  	//testDataSet.writePointData(t, mpiComm);
	  
	  	testDataSet.writeTimestepAttributes();
	  	testDataSet.writeVariables();

	  
	  	testDataSet.closeFile();

	  	// Delete 
		if (_x_data != NULL)
			delete [] _x_data;
		if (_y_data != NULL)
			delete [] _y_data;
		if (_z_data != NULL)
			delete [] _z_data;
		if (_pressure_data != NULL)
			delete [] _pressure_data;


		logStream << "\nData creation for timestep    " << t << " took ($): " << dataCreation.getDuration() << " s. " << std::endl;
		logStream << "h5hut conversion for timestep " << t << " took (^): " << dataConversion.getDuration() << " s. " << std::endl << std::endl << std::endl;

		// Sync writing before going to the next step
		MPI_Barrier(mpiComm);
	}
}



inline void ScalingTest::runWriteTest(int timeStep, size_t _numParticles)
{
	Timer dataWriting;

	for (int iter=0; iter<iterationCount; iter++)
	{
	  dataWriting.start();
		//testDataSet.writePointData(timeStep, mpiComm);
	  dataWriting.stop();

	  	logStream << iter << " ~ h5hut writing "<< _numParticles << " particles, timestep: " << timeStep << " took(#) " << dataWriting.getDuration() << " s."<< std::endl;
	}
}




#endif