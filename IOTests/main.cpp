#include <iostream>

#include "hdf5SimIO.h"

int main(int argc, char * argv[])
{
	Flecsi_Sim_IO::HDF5SimIO testDataSet( argv[1] );


	double *pressureData = new double[4*2*2];

	pressureData[0] = 10.0;		pressureData[1] = 20.0;		pressureData[2] = 30.0; 	pressureData[3] = 40.0;
	pressureData[4] = 12.0;		pressureData[5] = 14.0;		pressureData[5] = 16.0; 	pressureData[7] = 18.0;
	pressureData[8] = 22.0;		pressureData[9] = 24.0;		pressureData[10] = 26.0; 	pressureData[11] = 28.0;
	pressureData[12] = 32.0;	pressureData[13] = 34.0;	pressureData[14] = 36.0; 	pressureData[15] = 38.0;


	Flecsi_Sim_IO::Variable testScalar;
	testScalar.name = "pressure";
	testScalar.numDims = 3;
	testScalar.gridDims.push_back(4);
	testScalar.gridDims.push_back(2);
	testScalar.gridDims.push_back(2);
	testScalar.varType = Flecsi_Sim_IO::grid_cellCentered;
	testScalar.dataType = "double";
	testScalar.data = pressureData;

	Flecsi_Sim_IO::TimeStep testTs;
	testTs.timestep = 0;
	testTs.timeStamp = 0.125;
	testTs.vars.push_back(testScalar);

	testDataSet.timeVariables.push_back(testTs);
	testDataSet.createDataset();
	testDataSet.writeGridData(0);

	return 0;
}
