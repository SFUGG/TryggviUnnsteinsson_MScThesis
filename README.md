# Modelling glaciovolcanic caves and chimneys - MSc Thesis (Simon Fraser University)
## Tryggvi Unnsteinsson

This folder contains the code we used to investigate the formation and melt-through of glaciovolcanic voids in chapter 3 of my MSc thesis.

`define_variables.sh` - A file that takes in the user inputs of the desired physical parameters, imports them to the appropriate files and preps the working directory for simulation.

`run_iterative.sh` - The file that runs the simulation and iterates between different parameter combinations. It is created from the file `run_iterative_fillin.sh` which is edited by `define_variables.sh` to contain the appropriate variables.

`elmer_job.sh` - Equivalent as previous file but is used to run the code on Compute Canada's supercomputer clusters. (This file will be edited to not contain account details before publishing.)

`create_geometry.py` - Creates the mesh for each parameter combination and void geometry. The mesh is created using [GMSH](https://gmsh.info/) and exported in `.msh` format.

`diangostic.sif` - The diagnostic solver input file for Elmer/Ice. The file is exported into each parameter folder with the appropriate parameters changed.

`prognostic.sif` - The prognostic solver input file for Elmer/Ice. The file is exported into each parameter folder with the appropriate parameters changed.

`remesh_deformed.py` - A file that takes in the prognostic `.vtu` result, remeshes the void and top surface and returns a new mesh in `.msh` format.

## Execution

To use the code, navigate to a working directory where the results are to be saved and download the "main" folder to said directory.

To setup for a simulation for desired physical parameters (glacier thickness, bed slope, total heat flux and mode of heat transfer) update the `define_variables.sh` bash file to include those parameters (see [Physical parameters](https://github.com/TryggviU/ElmerIce_GCC/blob/main/README.md#physical-parameters)). Each variable has to have at least one option specified. Once the variables have been defined run the command:
```
bash main/define_variables.sh
```
This should create parameter folders with names of the format 'H{}\_I{}\_Q{}\_M{}', e.g. for 100 m thickness, 5° bed slope, 2500 kW total heat flux and mode of heat flux nr. 3 the associated folder name is 'H100\_I5\_Q2500\_M3'. A `.txt` file containing all the different parameter folder names is also created. The command should also create two separate files within the `main` directory, `elmer_job.sh` and `run_iterative.sh`, which are used to run the code on either a computer cluster or your personal computer, respectively.

Then to run the simulation over the chosen parameters simply run the command:
```
bash main/run_iterative.sh 
```
Simulations with each parameter composition should then run in serial and populate the parameter folders with the `.vtu` result files from each iteration which can be viewed in [Paraview](https://www.paraview.org/).

## Physical parameters

The user has four variables to choose from:
1. Glacier thickness, *H* [m].
   - The thicknesses used in our simulations are 50, 100, 150 and 200 meters.
3. Bed inclination, *I* [°].
   - The bed slopes used in our simulations are 0°, 5°, 10° and 15°.
4. Total heat flux, *Q* [kW].
   - The heat fluxes used in our simulations range from 0 to 3.5 kW.
5. The mode of heat transfer, *M*.
   1. Uniform melt over the entire surface.
   2. Melt is a function of radial distance from vertical.
   3. Melt is a function of radial distance from the centre line of the void.
   4. Melt is a function of radial distance from the geothermal point source.
   
   Modes (ii) and (iii) compute a cylindroid surface with some minimal radial distance (R_min) that the total heat flux is applied to. Specific heat fluxes decrease radially, q_M = f(R_min/R), from the surface of that cylindroid. If a point is within the cylindroid the specific heat flux is the same as at the surface. The melt rate (dr/dt) is computed as: dr/dt = q / (rho_i * L_f). The melt amount is hence: dr = q / (rho_i * L_f) * dt.
  
Additional variables which the user must define are:
- Geometry (geo), as either "cave" or "chimney".
- Initial geometry radius (R).
- Plume radius (r)
- Timestep of the prognostic iterations in hours (dt).
- Total number of iterations (N).
- Total runtime for each simulation in hours (hrs).
- Maximum number of simultaneous simulations if run on a cluster (max_number_simulations)
