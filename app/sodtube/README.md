![logo](doc/flecsph_logo_bg.png)

# Sod Shock Tube test

This test implements five standard Sod shock tube tests in 1D 
(see e.g. Toro 1999, "Riemann Solvers and Numerical Methods for Fluid Dynamics"). 

## Generating initial data 
Use the generator as follows: 

    % mpirun -np X ./sodtube_generator [-n <number of particles> [-t <testnum>]]

where `testnum` is a number from 1 to 5 (1 by default), and number of particles
is 1000 by default. 
This will produce an h5part file named `h5part_sodtube.h5part` to be input into 
the evolution code. 

## Running the evoluiton app 

    % mpirun -np X ./sodtube 

As long as FleCSI does not provide a way to read the file name, it is hardcoded 
in the main_driver.cc file to be `h5part_sodtube.h5part`. 

The current version generates a file `output_sodtube.h5part` with the result. 

## Visualize resutls

Paraview can be used after loading the h5part module. 
