#include <iostream>
#include <string>

#include "H5Cpp.h"



const H5std_string	FILE_NAME("h5tutr_dset.h5");
const H5std_string	DATASET_NAME("dset");
const int	 NX = 4;                     // dataset dimensions
const int	 NY = 6;
const int	 NZ = 4;
const int	 RANK = 2;

int main (int argc, char *argv[])
{
	if (argc <2)
	{
		std::cout << "Filename needed!. Exiting now ..." << std::endl;
		return 0;
	}


    try
    {
		// Turn off the auto-printing when failure occurs so that we can
		// handle the errors appropriately
		H5::Exception::dontPrint();
		hid_t file, status;

		
		// Create a new file using the default property lists. 
		file = H5Fcreate(argv[1], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		
		hid_t    dataset, datatype, dataspace;
		// Create the data space for the dataset.
		hsize_t dims[3];              
		//hsize_t dims[2];               
		dims[0] = NX;
		dims[1] = NY;
		dims[2] = NZ;
		dataspace = H5Screate_simple(RANK, dims, NULL);


		
		//Define a datatype for the data in the dataset.
		//We will store little endian integers.
		
		datatype = H5Tcopy(H5T_NATIVE_INT);
		status = H5Tset_order(datatype, H5T_ORDER_LE);
		
		// Create a new dataset within the file using the defined
		// dataspace and datatype and default dataset creation
		// properties.
		// NOTE: H5T_NATIVE_INT can be used as the datatype if
		// conversion to little endian is not needed.
		
		dataset = H5Dcreate(file, "dataset_test", datatype, dataspace,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		H5Tclose(datatype);
		H5Dclose(dataset);
		H5Fclose(file);
    }
    catch(H5::FileIException error)			// catch failure caused by the H5File operations
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)		// catch failure caused by the DataSet operations
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)	// catch failure caused by the DataSpace operations
    {
		error.printError();
		return -1;
    }

    return 0;