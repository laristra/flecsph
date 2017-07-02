![logo](doc/flecsph_logo_bg.png)

# Sod Shock Tube test

This test implement the 1D sod shock tube test. 

## Generate the data 
Using the generator like: 

    % mpirun -np X ./sodtube_generator <number of particles>

This will generates a h5part file name h5part_sodtube.h5part which is directly
read by the program. 

## Running the application 

    % mpirun -np X ./sotube 

As long as FleCSI does not provide a way to read the file name, it is hardcoded 
in the main_driver.cc file. 

The current version generate an output_sodtube.h5part file with the result. 

## Visualize resutls

Paraview can be used after loading the h5part module. 
