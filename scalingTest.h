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

  public:
  	ScalingTest(){ theadingON=false; myRank=0; numRanks=1; iterationCount=1; }
  	~ScalingTest(){};

  	void enableTheading(){ theadingON=true; }
  	void setIterationCount(int _iterationCount){ iterationCount=_iterationCount; }
  	void initMPIScaling(MPI_Comm _comm);

  	void runScalingTest(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename);
  	void createPseudoData(size_t _numParticles, int numTimesteps, Octree simOct, std::string filename);
  	void runWriteTest(int timeStep, size_t _numParticles);
  	void readDatasetTest(std::string filename);

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

	  	testDataSet.setTimeStep(t);


	  	// Add timestep attributes
	  	Flecsi_Sim_IO::Attribute timeValue("sim_time", Flecsi_Sim_IO::timestep, "float", (float)t/10.0);
	  	testDataSet.timestepAttributes.push_back(timeValue);

	  	testDataSet.addTimeStepAttribute( Flecsi_Sim_IO::Attribute("step_time", Flecsi_Sim_IO::timestep, "int32_t", t) );

	  	

	  	//
	  	// Add Variables

	  	// Option 1: create and push to timstep
	  	Flecsi_Sim_IO::Variable _x("x", Flecsi_Sim_IO::point, "float", myNumParticles, _x_data);
	  	testDataSet.vars.push_back(_x);

	  	// Option 2: init, create and push
		Flecsi_Sim_IO::Variable _y;
		_y.createVariable("y", Flecsi_Sim_IO::point, "float", myNumParticles, _y_data);
		testDataSet.vars.push_back(_y);

		// Option 3: create and add in 2 steps
		Flecsi_Sim_IO::Variable _z("z", Flecsi_Sim_IO::point, "float", myNumParticles, _z_data);
		testDataSet.addVariable( _z );

		// Option 4: create and add, all in one step
		testDataSet.addVariable( Flecsi_Sim_IO::Variable("pressure", Flecsi_Sim_IO::point, "double", myNumParticles, _pressure_data) );

	  dataConversionTimer.stop();
	  

	  	// Write data to hard disk
	  dataWritingTimer.start();
	  	testDataSet.writeTimestepAttributes();
	  	testDataSet.writeVariables();
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


inline void ScalingTest::readDatasetTest(std::string filename)
{
	Flecsi_Sim_IO::HDF5ParticleIO testDataSet( filename.c_str(), Flecsi_Sim_IO::READING, mpiComm );


	std::cout << "# dataset attributes:" << testDataSet.getNumDatasetAttributes() << ", " << testDataSet.datasetAttributes.size() << std::endl;

	for (int i=0; i<testDataSet.getNumDatasetAttributes(); i++)
	{
		std::cout << testDataSet.datasetAttributes[i].name;
		std::string type = testDataSet.datasetAttributes[i].dataType;
		if (type == "int32_t")
		{
			std::cout << " : " <<testDataSet. datasetAttributes[i].getAttributeValue<int32_t>() << std::endl;
		}
		else if (type == "int64_t")
		{
			std::cout << " : " << testDataSet.datasetAttributes[i].getAttributeValue<int64_t>() << std::endl;
		}
		else if (type == "float")
		{
			std::cout << " : " << testDataSet.datasetAttributes[i].getAttributeValue<float>() << std::endl;
		}
		else if (type == "double")
		{
			std::cout << " : " << testDataSet.datasetAttributes[i].getAttributeValue<double>() << std::endl;
		}
		else if (type == "string")
		{
			std::cout << " ::: "  << std::endl;
			char * _temp = new char[ testDataSet.datasetAttributes[i].numElements ];
			testDataSet.datasetAttributes[i].getAttributeArray<char>(_temp);
			std::cout << " : " << _temp << std::endl;
		}
	}
}


/*
inline void ScalingTest::readDatasetTest(std::string filename)
{
	// 
	Flecsi_Sim_IO::HDF5ParticleIO testDataSet( filename.c_str(), Flecsi_Sim_IO::READING, mpiComm );


	// Get information about file
	//if (myRank == numRanks-1)
  

	if (myRank == 0)
	{
		// Dataset Attributes
		std::cout << myRank << " ~ #attributes: " << testDataSet.getNumDatasetAttributes() << std::endl;
		for (int i=0; i<testDataSet.getNumDatasetAttributes(); i++)
		{
			std::string varname, type;
			int numElements;
			testDataSet.getDatasetAttributeInfo(i, varname, type, numElements);
			

			std::cout << i << " : " << varname << ", " << type << ", " << numElements << "   ";

			// Using templates
			if (i == 0)
			{
				int64_t _attrib;
				_attrib = testDataSet.readDatasetAttribute<int64_t>(varname, type);
				std::cout << "value: " << _attrib << std::endl;
			}
			else
			{
				// using void* pointers
				void *_attrib;
				if (type == "int32_t")
				{
					_attrib = new int32_t[numElements];
					testDataSet.readDatasetAttributeArray(varname, type, (int32_t *)_attrib);
					std::cout << "value: " << ((int32_t *)_attrib)[0] << std::endl;
				}
				else if (type == "int64_t")
				{
					_attrib = new int64_t[numElements];
					testDataSet.readDatasetAttributeArray(varname, type, (int64_t *)_attrib);
					std::cout << "value: " << ((int64_t *)_attrib)[0] << std::endl;
				}
				else if (type == "float")
				{
					_attrib = new float[numElements];
					testDataSet.readDatasetAttributeArray(varname, type, (float *)_attrib);
					std::cout << "value: " << ((float *)_attrib)[0] << std::endl;
				}
				else if (type == "double")
				{
					_attrib = new double[numElements];
					testDataSet.readDatasetAttributeArray(varname, type, (double *)_attrib);
					std::cout << "value: " << ((double *)_attrib)[0] << std::endl;
				}
				else if (type == "string")
				{
					_attrib = new char[numElements];
					testDataSet.readDatasetAttributeArray(varname, type, (char *)_attrib);
					std::cout << "value: " << ((char *)_attrib) << std::endl;
				}
			}
		}

		std::cout << myRank << " ~ #timesteps: "  << testDataSet.getNumTimesteps() << std::endl;
		
		for (int i=0; i<testDataSet.getNumTimesteps(); i++ )
		{
			std::cout << myRank << " ~ timestep: "    << i 
								<< " : #attributes: " << testDataSet.getNumTimestepAttributes(i) 
								<< ", #variables: "   << testDataSet.getNumVariables(i) 
								<< ", #particles: "   << testDataSet.getNumPartcles(i) << std::endl;

			for (int j=0; j<testDataSet.getNumVariables(i); j++)
			{
				std::string varname, type;
				int numElements;
				testDataSet.getTimestepVariableInfo(j, varname, type, numElements);
				std::cout << j << " : " << varname << ", " << type << ", " << numElements << std::endl;
			}
		}
	}
  
	MPI_Barrier(mpiComm);

	for (int i=0; i<testDataSet.getNumTimesteps(); i++ )
	{
		for (int j=0; j<testDataSet.getNumVariables(i); j++)
		{
			std::string varname, type;
			int numElements;
			testDataSet.getTimestepVariableInfo(j, varname, type, numElements);
			std::cout << j << " : " << varname << ", " << type << ", " << numElements << std::endl;

			if (i == 0 && j == 0)
			{
				//double *_data = new double[numElements/numRanks];
				double *_data = new double[numElements];
				testDataSet.readVariable(varname, type, _data);
				for (int k=0; k<numElements/numRanks; k++)
				{
					std::cout << myRank << " ~ " << k << ": " << ((double *)_data)[k] << std::endl;
				}
			}
		}
	}
}
*/

#endif