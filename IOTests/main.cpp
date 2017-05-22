#include <iostream>

#include "hdf5SimIO.h"

int main(int argc, char * argv[])
{

 	#if 1
	//
	// Point data test
	

	// Some variables
	float *_x_data = new float[16];
	float *_y_data = new float[16];
	float *_z_data = new float[16];
	float *_pressure_data = new float[16];

	_x_data[0] = 10.0; 	_y_data[0] = 10.0; 	_z_data[0] =  8.0; 		_pressure_data[0] = 0.7;
	_x_data[1] =  5.0; 	_y_data[1] = 20.0; 	_z_data[1] = 12.0; 		_pressure_data[1] = 0.1;
	_x_data[2] = 20.0; 	_y_data[2] = 30.0; 	_z_data[2] = 21.0; 		_pressure_data[2] = 0.5;
	_x_data[3] =  0.0; 	_y_data[3] =  1.0; 	_z_data[3] = 27.0; 		_pressure_data[3] = 0.9;
	_x_data[4] =  7.0; 	_y_data[4] =  7.0; 	_z_data[4] =  2.0; 		_pressure_data[4] = 0.4;
	_x_data[5] = 11.0; 	_y_data[5] =  8.0; 	_z_data[5] =  8.0; 		_pressure_data[5] = 0.2;
	_x_data[6] = 26.0; 	_y_data[6] = 25.0; 	_z_data[6] = 15.0; 		_pressure_data[6] = 0.6;
	_x_data[7] = 31.0; 	_y_data[7] = 17.0; 	_z_data[7] = 18.0; 		_pressure_data[7] = 0.5;
	_x_data[8] =  4.0; 	_y_data[8] = 31.0; 	_z_data[8] = 30.0; 		_pressure_data[8] = 0.2;
	_x_data[9] =  8.0; 	_y_data[9] = 12.0; 	_z_data[9] =  9.0; 		_pressure_data[9] = 0.3;
	_x_data[10] = 10.0; _y_data[10] = 29.0; _z_data[10] = 24.0; 	_pressure_data[10] = 0.9;
	_x_data[11] = 10.0; _y_data[11] = 31.0; _z_data[11] = 12.0; 	_pressure_data[11] = 0.8;
	_x_data[12] = 30.0; _y_data[12] = 30.0; _z_data[12] =  7.0; 	_pressure_data[12] = 0.8;
	_x_data[13] = 24.0; _y_data[13] = 12.0; _z_data[13] =  8.0; 	_pressure_data[13] = 0.1;
	_x_data[14] = 16.0; _y_data[14] = 16.0; _z_data[14] = 14.0; 	_pressure_data[14] = 0.2;
	_x_data[15] =  8.0; _y_data[15] = 18.0; _z_data[15] = 31.0; 	_pressure_data[15] = 0.3;


	Flecsi_Sim_IO::Variable _x;
	_x.name = "x";
	_x.varType = Flecsi_Sim_IO::point;
	_x.dataType = "float";
	_x.data = _x_data;

	Flecsi_Sim_IO::Variable _y;
	_y.name = "y";
	_y.varType = Flecsi_Sim_IO::point;
	_y.dataType = "float";
	_y.data = _y_data;

	Flecsi_Sim_IO::Variable _z;
	_z.name = "z";
	_z.varType = Flecsi_Sim_IO::point;
	_z.dataType = "float";
	_z.data = _z_data;

	Flecsi_Sim_IO::Variable _pressure;
	_pressure.name = "pressure";
	_pressure.varType = Flecsi_Sim_IO::point;
	_pressure.dataType = "double";
	_pressure.data = _pressure_data;

	//
	// Create hdf5
	Flecsi_Sim_IO::HDF5SimIO testDataSet( argv[1] );

	testDataSet.createDataset();

	testDataSet.numDims = 3;
	testDataSet.datasetType = Flecsi_Sim_IO::structuredGrid;

	testDataSet.vars.push_back(_x);
	testDataSet.vars.push_back(_y);
	testDataSet.vars.push_back(_z);
	testDataSet.vars.push_back(_pressure);

	testDataSet.writePointData();

  #else
	// 
	// Grid data test
	
	//
	// Create data
	double *pressureData = new double[4*2*2];

	pressureData[0] = 10.0;		pressureData[1] = 20.0;		pressureData[2] = 30.0; 	pressureData[3] = 40.0;
	pressureData[4] = 12.0;		pressureData[5] = 14.0;		pressureData[5] = 16.0; 	pressureData[7] = 18.0;
	pressureData[8] = 22.0;		pressureData[9] = 24.0;		pressureData[10] = 26.0; 	pressureData[11] = 28.0;
	pressureData[12] = 32.0;	pressureData[13] = 34.0;	pressureData[14] = 36.0; 	pressureData[15] = 38.0;


	Flecsi_Sim_IO::Variable testScalar;
	testScalar.name = "pressure";
	testScalar.varType = Flecsi_Sim_IO::grid_cellCentered;
	testScalar.dataType = "double";
	testScalar.data = pressureData;


	//
	// Create hdf5
	Flecsi_Sim_IO::HDF5SimIO testDataSet( argv[1] );

	testDataSet.createDataset();
	testDataSet.vars.push_back(testScalar);
	
	testDataSet.datasetType = structuredGrid;
	testDataSet.numDims = 3;
	testDataSet.gridDims.push_back(4);
	testDataSet.gridDims.push_back(2);
	testDataSet.gridDims.push_back(2);

	testDataSet.writeGridData();
  #endif

	return 0;
}


// Run:
// ./hdf5Test gridTest.h5