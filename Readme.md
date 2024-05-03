# N-Particle System with Pairwise Interaction on the 2D Torus.

## Summary
This repository holds the simulation of the interacting particle system (IPS) for the work
**insert paper title here**.  
It holds the C++ source code that can be inspected and easily modified, an executable (compiled on Linux) simulating a single trajectory of the IPS and printing results to files, and two simple Python scripts to visualize the results.

## Introduction
**brief summary of the system**

## Build From Source
To create an executable on your platform, download the files in the src folder and compile the main.cpp file. Due to the use of OpenMP, a corresponding flag will need to be passed, depending on the compiler. It might also be necessary to specify the C++17 standard (or higher). On Linux, using gcc, compilation is invoked via `g++ -fopenmp -O3 -std=c++17 -o IPS.out main.cpp`.  
If you see an error that the file `cxxopts.hpp` has not been found, download it from https://github.com/jarro2783/cxxopts/tree/master/include and store in your system's include folder (on Linux, this is typically /usr/include). 

## Quickstart   
To run a simulation, run the executable with `./IPS.out`. This will simulate a single trajectory of an IPS with default parameters. It will create a .csv output file holding time series data of the mean
distance to the center of mass (COM), the mean squared displacement (MSD), and the (instantaneous) kinetic temperature (Tkin).  
The properties of the simulation can be controlled via various flags (eg. stepsize, number of particles, number of iterations, which integrator to use etc.). Run `./IPS.out --help` to view the possible options.
  
To plot the results, run `python plot_data.py <file>`. It accepts a flag `--title plot_title` to specify a title for the plot. Run `python plot_data.py --help` for more info.
  
In order to create an animation of a trajectory, `./IPS.out` needs to be executed with the `--trajectory` flag which will prompt the program to print out trajectory data to a second file.  
This file can be read with the second Python script via `python plot_trajectory <file>`. Once again, there are various options available (view them via the `--help` flag).
  
A snapshot of an animation can be seen below.
.![snapshot_animation_github](https://github.com/SchroedingersLion/Interacting_Particle_Systems/assets/70909827/b1cfc6b6-a946-442a-bb4f-a3bdbe877f05)
## Add additional forcefields or observables
**explain how users can add new force fields or change the observables that are collected**
