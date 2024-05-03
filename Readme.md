# N-Particle System with Pairwise Interaction on the 2D Torus.

## Summary
This repository holds the simulation of the interacting particle system (IPS) for the work
**insert paper title here**.  
It holds the C++ source code that can be inspected and easily modified, an executable simulating a single trajectory of the IPS and printing results to files, and two simple Python scripts to visualize the results.

## Introduction
**brief summary of the system**

## Quickstart
To run a simulation, run the executable with `./IPS.exe`. This will simulate a single trajectory of an IPS with default parameters. It will create a .csv output file holding time series data of the mean
distance to the center of mass (COM), the mean squared displacement (MSD), and the (instantaneous) kinetic temperature (Tkin).  
The properties of the simulation can be controlled via various flags (eg. stepsize, number of particles, number of iterations, which integrator to use etc.).  
Run `./IPS.exe --help` to view the possible options.
  
To plot the results, run `python plot_data.py <file>`. It accepts a flag `--title plot_title` to specify a title for the plot. Run `python plot_data.py --help` for more info.
  
In order to create an animation of a trajectory, `./IPS.exe` needs to be executed with the `--trajectory` flag which will prompt the program to print out trajectory data to a second file.  
This file can be read with the second Python script via `python plot_trajectory <file>`. Once again, there are various options available (view them via the `--help` flag).
  
A snapshot of an animation can be seen below.

## Add additional forcefields or observables
**explain how users can add new force fields or change the observables that are collected**.