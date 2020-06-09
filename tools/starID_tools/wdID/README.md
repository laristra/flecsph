## Files to set up, relax, and evolve a White Dwarf with 

The simulations should be done with the latest version of FleCSPH on the master branch. 

### General Instruction

1.  Create a directory, e.g. `white_dwarf`, on scratch where you will run your simulations. Copy all data files (and python scripts) from this directory into `white_dwarf`
2.  Create a backup of `wd_radial_profile_rho01e7.dat` and remove the header lines and last line from the original file. You should only have a file with the data columns at the end
3.  Copy the following executables from the flecsph directory: `sedov_3d_generator`, `hydro_3d`, `newtonian_3d`
4.  Get an interactive allocation 
5.  Load all the modules required to run `flecsph`
4.  Initialize the white dwarf by running e.g. (this example is for one mpi process)
    ```{engine=sh}
        mpirun -np 1 ./sedov_3d_generator initial.par
    ```
5.  Check what file(s) were created. Visualize the resulting h5part file with Paraview. 
6.  You can also run the script `1D_slice.py` to get a quick view of the density vs. radial distance for each particle. For that, run 
    ```{engine=sh}
        ./1D_slice.py -ifile wd_initial.h5part 
    ```
    and plot the resulting `1d_output_00000.dat` with gnuplot, e.g. 
    ```{engine=sh}
        gnuplot '1d_output_00000.dat' using 4:5 with points pt 7 ps 0.8
    ```
7.  As you can see from the output (either in paraview or gnuplot) the particles are distributed on a regular spherical lattice. If we run the simulation with this initial data
    the particles will most likely rearrange themselves in a more
    natural and irregular distribution. To speed up this process, we
    run a relaxation step without gravity. Instead,
    we create an effective potential, based on the radial density profile `wd_radial_profile_rho01e7.dat` and let the particles evolve in it. Ideally this should be done until the 
    particle density distribution agrees with the initial profile and does not change anymore. For 125,000 particles this takes >300,000 iterations. We will have to see how 
    feasible this is for >1,000,000 particles. However, to run the relaxation stage, use the following (example is again for one mpi process):
    ```{engine=sh}
        mpirun -np 1 ./hydro_3d relaxation.par
    ```    
8.  You can restart the simulation if you run out of time during relaxation (see below)
9.  To see if the white dwarf reached a relaxed state, you can use paraview or (better) a tool for the density vs. radial distance, e.g. `1D_slice.py` or a python script in the 
    flecsph tools folder 
10. Once the white dwarf is relaxed, evolve it without the external potential but with self-gravity via FMM by running the following commands (example for one mpi process):
    ```{engine=sh}
        mpirun -np 1 ./newtonian_3d evolution.par 
    ```

### Create new white dwarf radial profiles as initial data:
    Run the `WDtov_solver.py` script via
    ```{engine=sh}
        ipython WDtov_solver.py arg1 arg2 arg3
    ```
    where 
    arg1: central density in g/cm^3, e.g. 1e7
    arg2: width of radial step, e.g. 1e6
    arg3: GR corrections: True or False


Restart during e.g. relaxation (or evolution) phase

 * Check the file `scalar_reductions.dat` for the last iteration output, e.g. 346600
 * Back up your original relaxation or evolution file as e.g. `wd_relaxation_it346600.h5part`
 * Use the last iteration number as the `initial_iteration` parameter in your parameter file
 * Use the name of your original relaxation or evolution h5part file in `initial_data_prefix`
 * Run your simulation as before, e.g. 
  ```{engine=sh}
      mpirun -np 1 ./hydro_3d relaxation.par
  ```
Below is some typical values for WD

```
Mtot: 1.59517068564e+33 g; Rtot 704000000.0 cm
rho_central: 10000000.0 g/cm^3
P_central: 8.37239119678e+23 g/cm^3
```
