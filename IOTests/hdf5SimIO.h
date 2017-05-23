#ifndef _HDF5_SIM_IO_H
#define _HDF5_SIM_IO_H

#include "simIO.h"
#include "H5Cpp.h"

#include <iostream>

namespace Flecsi_Sim_IO
{

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

class HDF5SimIO: public SimIO
{
  public:
    HDF5SimIO():SimIO(){}
    HDF5SimIO(std::string _outputFileName):SimIO(_outputFileName){}

    hid_t getHDF5Datatype(std::string _datatypeName);

    int createDataset();
    int writeGridData();
    int writePointData();
};



hid_t HDF5SimIO::getHDF5Datatype(std::string _datatypeName)
{
    hid_t _datatype;
    if (_datatypeName == "double")
        _datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    else if (_datatypeName == "int")
        _datatype = H5Tcopy(H5T_NATIVE_INT);
    else if (_datatypeName == "float")
        _datatype = H5Tcopy(H5T_NATIVE_FLOAT);
    else if (_datatypeName == "char")
        _datatype = H5Tcopy(H5T_NATIVE_CHAR);
    else
    {
        std::cout << "Datatype " << _datatypeName << " has not been defined yet!" << std::endl;
        _datatype = NULL;
    }

    return _datatype;
}


inline int HDF5SimIO::createDataset()
{
	// Create file, will fail if already exists
	hid_t dataFile = H5Fcreate(outputFileName.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	H5Fclose(dataFile);

	return 0;
}



inline int HDF5SimIO::writeGridData()
{
	hid_t  dataFile;    
    dataFile = H5Fopen(outputFileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    for (int i=0; i<vars.size(); i++)
    {
    	if ( vars[i].varType == grid_cellCentered )
    	{
            hid_t  dataset, datatype;
            herr_t status;

            std::string _varname  = vars[i].name;

    		int _numDims  = numDims;
    		hsize_t *dims = new hsize_t[_numDims];
            for (int j=0; j<_numDims; j++)
                dims[j]=gridDims[j];

    		datatype = getHDF5Datatype( vars[i].dataType );
    		status   = H5Tset_order(datatype, endian==little?H5T_ORDER_LE:H5T_ORDER_BE);


            hid_t dataspace = H5Screate_simple(_numDims, dims, NULL);
    		dataset = H5Dcreate(dataFile, _varname.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    		status  = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, vars[i].data);

            //
            // Clean up
    		delete []dims;

            H5Sclose(dataspace);
            H5Tclose(datatype);
            H5Dclose(dataset);
    	}	
    }

    H5Fclose(dataFile);
	return 0;
}


inline int HDF5SimIO::writePointData()
{
    hid_t  dataFile;    
    dataFile = H5Fopen(outputFileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    for (int i=0; i<vars.size(); i++)
    {
        if ( vars[i].varType == point )
        {
            hid_t  dataset, datatype;
            herr_t status;

            std::string _varname  = vars[i].name;
            int _numDims          = numDims;

            datatype = getHDF5Datatype( vars[i].dataType );
            status   = H5Tset_order(datatype, endian==little?H5T_ORDER_LE:H5T_ORDER_BE);


            hid_t dataspace = H5Screate(H5S_SCALAR);
            dataset = H5Dcreate(dataFile, _varname.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status  = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, vars[i].data);

            H5Sclose(dataspace);
            H5Tclose(datatype);
            H5Dclose(dataset);
        }   
    }

    H5Fclose(dataFile);
    
    return 0;
}


} // end Flecsi_Sim_IO namespace

#endif