Libraries under consideration:
 1. HDF5
 2. H5hut
 3. Parallel-netcdf
 4. netcdf

H5Part - superceded by H5hut


Storage format: HDF5
Why? Lots of support available for vis + db + python scripts


What needs to be done:
 - Decide on Writing/Reading Methodology
 - Specify file format for output + checkpointing


Built:
 - hdf5				: /home/pascal/NGC/IO/hdf5-1.8.18/install
 - phdf5			: /home/pascal/NGC/IO/phdf5-1.8.18/install
 - H5hut			: /home/pascal/NGC/IO/H5hut-1.99.13/install
 - Parallel-netcdf	: /home/pascal/NGC/IO/parallel-netcdf-1.8.1/install


Note:
 - H5hut
 	- commented the internal of fn hdf5_set_fapl_mpiposix_property @line 929 in h5_hdf5_private.h and built on top of phdf5
	- coz hdf5 does not support posix anymore



