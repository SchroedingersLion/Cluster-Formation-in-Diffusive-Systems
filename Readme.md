# SimIPS.

## Summary
This repository holds a multithreaded C++ code for the **Sim**ulation of **I**nteracting **P**article **S**ystems (SimIPS) for the work
**Emergence of Clusters in Langevin Dynamics**.  
It holds the source code that can be easily modified, compiled and run, and two simple Python scripts to visualize the results.

## Introduction
We study a system of $N$ particles, whose configuration is given by the particle positions $\boldsymbol{x}_i\in \mathbb{R}^d$ and velocities $\boldsymbol{v}_i \in \mathbb{R}^d$, for $i=1,...,N$ and $d\in \{1,2,3\}$. The particles evolve in time $t$ via underdamped Langevin Dynamics, given by
$$
\begin{aligned}
\text{d}\boldsymbol{x}_i &=  \boldsymbol{v}_i \text{d}t, \\
\text{d}\boldsymbol{v}_i &= -\frac{1}{N} \sum_{i=j}^{N}\nabla W(\boldsymbol{x}_i, \boldsymbol{x}_j)\text{d}t -  \gamma \boldsymbol{v}_i \text{d}t + \sqrt{2\gamma \beta^{-1}}\text{d}\boldsymbol{B}_t^i,
\end{aligned}
$$
with friction $\gamma \ge 0$, temperature $\beta^{-1} \ge 0$ and a $d$-dimensional standard Wiener process $\boldsymbol{B}^i_t$.
Each particle pair interacts via interaction potential $W(\boldsymbol{x}, \boldsymbol{y})$, which is purely attractive in our case. For a list of implemented potentials, see below.  

The equations have to be solved numerically through suitable integrators, two of which are currently implemented (see the paper for a discussion).

For temperatures smaller than some critical temperature, i.e., $\beta$ larger than some critical $\beta_c$, the equilibrium distribution of particles is given by a single particle cluster which holds all particles. Starting the simulation with the initial positions sampled from a uniform distribution, for $\beta > \beta_c$ one observes one or more clusters forming. In the case of multiple clusters, the clusters randomly move around until they are close enough to another cluster such that the attractive forces leads to cluster merging, successively reducing the total number of clusters until the final equilibrium cluster remains.


## Build From Source
To create an executable on your platform, download the files in the src folder and compile the `main.cpp` file. Due to the use of OpenMP, a corresponding flag will need to be passed, depending on the compiler. It might also be necessary to specify the C++17 standard (or higher). On Linux, using gcc, compilation is invoked via `g++ -fopenmp -O3 -std=c++17 -o simips main.cpp`.  
If you see an error that the file `cxxopts.hpp` has not been found, download it from https://github.com/jarro2783/cxxopts/tree/master/include and store in your system's include folder (on Linux, this is typically /usr/include). 

## Quickstart   
To run a simulation, run the executable with `./simips`. This will simulate a single trajectory of an IPS with default parameters. It will create a `.csv` output file holding time series data of the mean distance to the center of mass (COM), the mean squared displacement (MSD), and the (instantaneous) kinetic temperature (Tkin).  
The properties of the simulation can be controlled via various flags (e.g., stepsize, number of particles, number of iterations, which integrator to use etc.). Run `./simips --help` to view the possible options.

The two Python scripts to visualize the results needs to be run from a Python environment (e.g., conda environment) with installed `numpy` and `matplotlib` packages.
To plot the time series results, run `python plot_data.py <file>`. It accepts a flag `--title plot_title` to specify a title for the plot. Run `python plot_data.py --help` for more info. 
  
In order to create an animation of a trajectory, `./simips` needs to be executed with the `--trajectory` flag, which will prompt the program to print out trajectory data to a second file.  
This file can be read with the second Python script via `python plot_trajectory <file>`. Once again, there are various options available (view them via the `--help` flag).
  
A snapshot of an animation can be seen below.
.![snapshot_animation_github](https://github.com/SchroedingersLion/Interacting_Particle_Systems/assets/70909827/b1cfc6b6-a946-442a-bb4f-a3bdbe877f05)

## Add additional forcefields or observables
**explain how users can add new force fields or change the observables that are collected**
