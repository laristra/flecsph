![logo](doc/flecsph_logo_bg.png)

# Kelvin-Helmholtz instabilities

The test implement Kelvin-Helmholtz instabilities. 
(see Modelling discontinuities and Kelvin-Helmholtz instabilities in SPH, 
Daniel J. Price)

## Generating initial data 
Use the generator as follows: 

    % mpirun -np X ./KH_generator ./parameter_file.par

## Running the evoluiton app 

  % mpirun -np X ./hydro_2d ./parameter_file.par

## Visualize resutls

Paraview can be used after loading the h5part module. 
