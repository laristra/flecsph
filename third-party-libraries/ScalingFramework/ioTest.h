#ifndef _IO_TESTS_H_
#define _IO_TESTS_H_


#include <string>
#include <iostream>
#include <sstream>
#include <random>

#include <mpi.h>

#include "timer.h"
#include "log.h"
#include "hdf5ParticleIO.h"
#include "simIO.h"


class IOTest
{
	int myRank;
	int numRanks;
	MPI_Comm mpiComm;
	std::stringstream logStream;

  public:
  	IOTest(){ myRank=0; numRanks=1; }
  	~IOTest(){};

  	void initMPI(MPI_Comm _comm);

  	void writeDatasetTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename);
  	void readDatasetTest(std::string filename);

  	std::string getTimingLog(){ return logStream.str(); }
};



inline void IOTest::initMPI(MPI_Comm _comm)
{
	mpiComm = _comm;
	MPI_Comm_rank(_comm, &myRank);
	MPI_Comm_size(_comm, &numRanks);
}



inline void IOTest::writeDatasetTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename)
{
	Timer dataCreationTimer, dataConversionTimer, dataWritingTimer;

	int myNumParticles = _numParticles/numRanks;

	//
	// Get Partition from octree
	std::string nodeID = simOct.getNodeAt( (size_t)myRank );

	float nodeExtents[6];
	simOct.getSpatialExtent(nodeID, nodeExtents);


	//
	// Create Dataset
	Flecsi_Sim_IO::HDF5ParticleIO testDataSet( filename.c_str(), Flecsi_Sim_IO::WRITING, mpiComm );



	testDataSet.writeDatasetAttribute("numParticles", "int64_t", 80);
	testDataSet.writeDatasetAttribute("gravitational constant", "double", 6.67300E-11);
	testDataSet.writeDatasetAttribute("numDims", "int32_t", 3);
	testDataSet.writeDatasetAttribute("use_fixed_timestep", "int32_t", 0);

	char *simName = "random_data";
	testDataSet.writeDatasetAttributeArray("sim_name", "string", simName);

	for (int t=0; t<numTimesteps; t++)
	{
		//
		// Generate fake timestep data
	  dataCreationTimer.start();
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
	  dataCreationTimer.stop();



	  	//
		// Push fake timestep data to h5hut
	  dataConversionTimer.start();

	  	Flecsi_Sim_IO::Timestep _ts;

	  	_ts.timeStep = t;


	  	// Add timestep attributes
	  	Flecsi_Sim_IO::Attribute timeValue("sim_time", "float", (float)t/10.0);
	  	_ts.addAttribute(timeValue);

	  	_ts.addAttribute( Flecsi_Sim_IO::Attribute("step_time", "int32_t", t) );


	  	//
	  	// Add Variables

	  	// Option 1: create and push to timstep
	  	Flecsi_Sim_IO::Variable _x("x", "float", myNumParticles, _x_data);
	  	_ts.addVariable(_x);

	  	// Option 2: init, create and push
		Flecsi_Sim_IO::Variable _y;
		_y.createVariable("y", "float", myNumParticles, _y_data);
		_ts.addVariable(_y);

		// Option 3: create and add in 2 steps
		Flecsi_Sim_IO::Variable _z("z", "float", myNumParticles, _z_data);
		_ts.addVariable(_z);

		// Option 4: create and add, all in one step
		_ts.addVariable( Flecsi_Sim_IO::Variable("pressure", "double", myNumParticles, _pressure_data) );

	  dataConversionTimer.stop();
	  

	  	// Write data to hard disk
	  dataWritingTimer.start();
	  	testDataSet.writeTimestep(_ts);
	  dataWritingTimer.stop();
	  

	  	// Clean up 
		if (_x_data != NULL)
			delete [] _x_data;
		if (_y_data != NULL)
			delete [] _y_data;
		if (_z_data != NULL)
			delete [] _z_data;
		if (_pressure_data != NULL)
			delete [] _pressure_data;

		logStream << "\nData creation for timestep    " << t << " took ($): " << dataConversionTimer.getDuration() << " s. " << std::endl;
		logStream << "h5hut conversion for timestep " << t << " took (^): " << dataConversionTimer.getDuration() << " s. " << std::endl << std::endl << std::endl;
		logStream << "h5hut writing "<< _numParticles << " particles, timestep: " << t << " took(#) " << dataWritingTimer.getDuration() << " s."<< std::endl;

		// Sync writing before going to the next step
		MPI_Barrier(mpiComm);
	}
}


inline void IOTest::readDatasetTest(std::string filename)
{
	Flecsi_Sim_IO::HDF5ParticleIO testDataSet( filename.c_str(), Flecsi_Sim_IO::READING, mpiComm );


	// Read database attributes
	if (myRank == 0)
	{
		std::cout << "# dataset attributes: " << testDataSet.getNumDatasetAttributes() << std::endl;

		for (int i=0; i<testDataSet.getNumDatasetAttributes(); i++)
		{
			std::cout << testDataSet.datasetAttributes[i].name;
			std::string type = testDataSet.datasetAttributes[i].dataType;
			if (type == "int32_t")
				std::cout << ": " <<testDataSet. datasetAttributes[i].getAttributeValue<int32_t>() << std::endl;
			else if (type == "int64_t")
				std::cout << ": " << testDataSet.datasetAttributes[i].getAttributeValue<int64_t>() << std::endl;
			else if (type == "float")
				std::cout << ": " << testDataSet.datasetAttributes[i].getAttributeValue<float>() << std::endl;
			else if (type == "double")
				std::cout << ": " << testDataSet.datasetAttributes[i].getAttributeValue<double>() << std::endl;
			else if (type == "string")
				std::cout << ": " << (char *)testDataSet.datasetAttributes[i].data << std::endl;
		}
		std::cout << "\n";
	}


	// Read timesteps
	int numTimesteps = testDataSet.getNumTimesteps();

	for (int t=0; t<numTimesteps; t++)
	{
		Flecsi_Sim_IO::Timestep _ts;

		testDataSet.readTimestepInfo(t, _ts);

		// display
		if (myRank == 0)
		{
			std::cout << "\n\nNum attributes: " << _ts.numAttributes << std::endl;
			for (int i=0; i<_ts.numAttributes; i++)
			{
				std::cout << _ts.attributes[i].name;
				std::string type = _ts.attributes[i].dataType;
				if (type == "int32_t")
					std::cout << ": " << _ts.attributes[i].getAttributeValue<int32_t>() << std::endl;
				else if (type == "int64_t")
					std::cout << ": " << _ts.attributes[i].getAttributeValue<int64_t>() << std::endl;
				else if (type == "float")
					std::cout << ": " << _ts.attributes[i].getAttributeValue<float>() << std::endl;
				else if (type == "double")
					std::cout << ": " << _ts.attributes[i].getAttributeValue<double>() << std::endl;
				else if (type == "string")
					std::cout << ": " << (char *)_ts.attributes[i].data << std::endl;
			}

			std::cout << "\nNum variables: " << _ts.numVariables << std::endl;
			for (int i=0; i<_ts.numVariables; i++)
				std::cout << i << " : " << _ts.vars[i].name << ", " << _ts.vars[i].dataType << ", " << _ts.vars[i].numElements << std::endl;
		}


		// Read data
		for (int i=0; i<_ts.numVariables; i++)
		{
			if (_ts.vars[i].dataType == "int32_t")
			{
				_ts.vars[i].data = new int32_t[_ts.vars[i].numElements];
				testDataSet.readVariable(_ts.vars[i].name, _ts.vars[i].dataType, (int32_t *)_ts.vars[i].data);
			}
			else if (_ts.vars[i].dataType == "int64_t")
			{
				_ts.vars[i].data = new int64_t[_ts.vars[i].numElements];
				testDataSet.readVariable(_ts.vars[i].name, _ts.vars[i].dataType, (int64_t *)_ts.vars[i].data);
			}
			else if (_ts.vars[i].dataType == "float")
			{
				_ts.vars[i].data = new float[_ts.vars[i].numElements];
				testDataSet.readVariable(_ts.vars[i].name, _ts.vars[i].dataType, (float *)_ts.vars[i].data);


				if (myRank == 0)
				{
					std::cout << "\n\n" << _ts.vars[i].name << ":\n";
					for (int k=0; k<_ts.vars[i].numElements/numRanks; k++)
						std::cout << myRank << " ~ " << k << ": " << ((float *)_ts.vars[i].data)[k] << std::endl;
				}
			}
			else if (_ts.vars[i].dataType == "double")
			{
				_ts.vars[i].data = new double[_ts.vars[i].numElements];
				testDataSet.readVariable(_ts.vars[i].name, _ts.vars[i].dataType, (double *)_ts.vars[i].data);

				if (myRank == 0)
				{
					std::cout << "\n\n" << _ts.vars[i].name << ":\n";
					for (int k=0; k<_ts.vars[i].numElements/numRanks; k++)
						std::cout << myRank << " ~ " << k << ": " << ((double *)_ts.vars[i].data)[k] << std::endl;
				}
			}
		}
	}
}


#endif