#ifndef _HDF5_SIM_IO_H
#define _HDF5_SIM_IO_H

#include "simIO.h"
#include "H5Cpp.h"

class HDF5SimIO: public SimIO
{
  public:
  	int createDataset();
	int writeData();
};



// 
// HDF5 Represenation

// Database: file
// Timestep: Group
// 		Each Variable: Dataset
//			Datatype: datatype
//			Type: Dataset attribute (point or grid or ...)
//				  Dataspace (for strunctured grid description with dims ...)
//			
//		Timestep : Group attribute
//		Timestamp: Group attribute


inline int createDataset()
{
	// Create file, will fail if already exists
	hid_t dataFile = H5Fcreate(outputFileName.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

	H5Fclose(dataFile);
	return 0;
}


inline int writeData(int _timeStep)
{
	hid_t  dataFile, fapl_id;        


    dataFile = H5Fopen(outputFileName.c_str(), H5F_ACC_RDWR, fapl_id);

    
    for (int i=0; i<timeVariables[_timeStep].size(); i++)
    {
    	hid_t  dataset, datatype, dataspace;
    	herr_t status;

    	if ( (timeVariables[_timeStep]).varType == point )
    	{
    		int _numDims = (timeVariables[_timeStep]).numDims;
    		std::string _varname = (timeVariables[_timeStep]).name;
    		
    		hsize_t *dims = new hsize_t[_numDims];


    		hid_t dataspace = H5Screate_simple(_numDims, dims, NULL);
    		datatype = H5Tcopy(H5T_NATIVE_INT);

    		status = H5Tset_order(datatype, endian==little?H5T_ORDER_LE:H5T_ORDER_BE);

    		dataset = H5Dcreate(dataFile, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    		status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    		delete []dims;
    	}

    	H5Sclose(dataspace);
    	H5Tclose(datatype);
    	H5Dclose(dataset);
    }

    H5Fclose(dataFile);
	return 0;
}

#endif