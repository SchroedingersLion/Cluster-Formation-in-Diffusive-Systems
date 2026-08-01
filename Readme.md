
# Cluster Formation in Diffusive Systems
This repository holds code to simulate diffusive N-particle systems. In particular, it holds the code for the numerical experiments in the work 
Leimkuhler, B., Lohmann, R., Pavliotis, G., Whalley, P.A., **Cluster formation for weakly interacting kinetic Langevin dynamics**. 2025. [https://arxiv.org/abs/2510.25034](https://arxiv.org/abs/2510.25034)

<img src="https://github.com/user-attachments/assets/d17d4966-8824-475a-afe3-50338c724f56" alt="example_snapshots_1D_2D" width="500">

It holds 
- src folder holding the C++ codebase **SimIPS** for the particle simulation. 
- Python script **plot_data.py** to plot the observable data printed from the simulation.
- Python script **plot_trajectory.py** to create an animation of the system trajectory printed from the simulation.
- Jupyter notebook **stability_and_fluctuation_analysis.ipynb** for the numerical treatment of the linear stability and fluctuation analysis in Sec. 3.
- Animations folder storing the example animations referred to in the paper.

The following sections briefly explain the usage of the code. Feel free to reach out with any questions via creating an issue.

## SimIPS
SimIPS (**Sim**ulation of **I**nteracting **P**article **S**ystems) is a lightweight C++ code run from the command line to simulate 
a system of $N$ particles in a cubic simulation box of edge length $L$ under periodic boundaries. It can treat one, two, or three dimensional systems. The particles interact via pairwise interaction potentials. The current implementation offers the three interaction potentials discussed in the paper. The system is then simulated by integrating Langevin dynamics (using one of two high-quality integrators), and observables are measured regularly. Their time series are printed to a csv file at the end. Optionally, the whole trajectory can be stored and printed to a file as well.  
The force calculation requires computing all pairwise distances between the particles, which is a $\mathcal{O}(N^2)$ operation---  
the computational bottleneck of the code. To speed this up, we use multithreading via OpenMP in the force routine. At the moment, the code does
not use advanced algorithmic approaches (such as cell or Verlet lists) to improve the scaling with $N,$ since these techniques would not have
helped much for the particular settings we studied. On a Dell Latitude 5530 with an i7-1265U processor, using 10 threads, running a system of 2,000 particles for 20,000 iterations leads to a runtime of 85 seconds.

It is possible to extend the codebase by adding observables to collect, pairwise interaction functions, and Langevin dynamics integrators.
To add or change the observables to collect, modify the "measurements.h" file. To add or change the pairwise interaction functions, modify the 
"model.h" file. To add new integrators, modify the "simulation.h" file. The files contain comments and instructions to explain more details.

### Build From Source
To create an executable on your platform, download the files in the "SimIPS/src" folder and compile the `main.cpp` file. Due to the use of OpenMP, a corresponding flag will need to be passed, depending on the compiler. It might also be necessary to specify the C++17 standard (or higher). On Linux, using gcc, compilation is invoked via  
`g++ -fopenmp -O3 -std=c++17 -o simips main.cpp`.  
This will create a `simips` executable.  
If you see an error that the file `cxxopts.hpp` has not been found, download it from  
[https://github.com/jarro2783/cxxopts/tree/master/include](https://github.com/jarro2783/cxxopts/tree/master/include) and store it in your system's include folder (on Linux, this is typically /usr/include). 

### Quickstart   
To run a simulation, run the executable with `./simips`. This will simulate a single trajectory of an IPS with default parameters. It will create a `results.csv` output file holding time series data of the mean distance to the centre of mass (COM), the mean squared displacement (MSD), and the (instantaneous) kinetic temperature (Tkin). The first column gives the corresponding simulation times. 
The properties of the simulation can be controlled via various flags (e.g., stepsize, number of particles, number of iterations, which integrator to use etc.). Run `./simips --help` to view the possible options.

The two Python scripts to visualise the results need to be run from a Python environment (e.g., conda environment) with installed `numpy` and `matplotlib` packages.
To plot the time series results, run `python plot_data.py <file>`. It accepts a flag `--title plot_title` to specify a title for the plot. Run `python plot_data.py --help` for more info.

In order to create an animation of a trajectory, `./simips` needs to be executed with the `--trajectory` flag, which will prompt the program to print out trajectory data to a second file.  
This file can be read with the second Python script via `python plot_trajectory <file>`. Once again, there are options available (view them via the `--help` flag).
  



