#include <iostream>

#include "hdf5SimIO.h"

int main(int argc, char * argv[])
{

	std::cout << "hi" << std::endl;
	Flecsi_Sim_IO::HDF5SimIO hdf5Output;

	hdf5Output testDataSet("testHDF5.h5");

	double *pressureData = new double[4*2*2];

	pressureData[0] = 10.0;		pressureData[1] = 20.0;		pressureData[2] = 30.0; 	pressureData[3] = 40.0;
	pressureData[4] = 12.0;		pressureData[5] = 14.0;		pressureData[5] = 16.0; 	pressureData[7] = 18.0;
	pressureData[8] = 22.0;		pressureData[9] = 24.0;		pressureData[10] =26.0; 	pressureData[11] = 28.0;
	pressureData[12] = 32.0;	pressureData[13] = 34.0;	pressureData[14] = 36.0; 	pressureData[15] = 38.0;


	Flecsi_Sim_IO::Variable testScalar;
	test.name = "pressure";
	test.numDims = 3;
	test.gridDims.push_back(4);
	test.gridDims.push_back(4);
	test.gridDims.push_back(2);
	test.varType = Flecsi_Sim_IO::grid_cellCentered;
	test.dataType = "double";
	test.data = pressureData;

	Flecsi_Sim_IO::TimeStep testTs;
	testTs.timestep = 0;
	testTs.timeStamp = 0.125;
	testTs.vars.push_back(testScalar);

	hdf5Output.timeVariables.push_back(testTs);
	hdf5Output.createDataset();
	hdf5Output.writeData();

	return 0;
}
