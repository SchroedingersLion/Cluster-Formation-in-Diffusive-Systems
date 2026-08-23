
# SimIPS: A C++ Simulation Engine for Interacting Particle Systems
[![DOI](https://zenodo.org/badge/750927167.svg)](https://zenodo.org/badge/latestdoi/750927167)

This repository contains `SimIPS`, a lightweight C++ engine for simulating pairwise interacting N-particle systems in a periodic simulation box. 

<img src="https://github.com/user-attachments/assets/d17d4966-8824-475a-afe3-50338c724f56" alt="example_snapshots_1D_2D" width="550">
<br>

<sub><b>Built with SimIPS:</b> Visualized cluster formation in a weakly interacting kinetic Langevin system.<br>
For the mathematical analysis of this type of cluster formation, see:<br>
<i>Leimkuhler, B., Lohmann, R., Pavliotis, G. A., & Whalley, P. A. (2026). Cluster formation for weakly interacting kinetic Langevin dynamics.</i> <a href="https://arxiv.org/abs/2510.25034">Arxiv</a><br>
📂 See also the <code>papers/Leimkuhler_et_al_2026/</code> directory.</sub>

## Installation

The repository contains the following structure:
- <code>src</code> folder holding the C++ codebase **SimIPS** for the particle simulation. 
- <code>scripts</code> folder holding Python scripts <code>plot_data.py</code> to plot the observable time series printed by **SimIPS** and <code>plot_trajectory.py</code> to create an animation of the system trajectory.
- <code>papers</code> folder holding the academic projects **SimIPS** was used in.

The following sections briefly explain the usage of the code. Feel free to reach out with any questions via creating an issue.

## SimIPS
SimIPS (**Sim**ulation of **I**nteracting **P**article **S**ystems) is a lightweight C++ code run from the command line to simulate 
a system of $N$ particles in a cubic simulation box of edge length $L$ under periodic boundaries. It can treat one, two, or three dimensional systems. The particles interact via pairwise interaction potentials $U(r_{i,j})$ with $r_{i,j}$ the Euclidean distance between particles $i$ and $j$. The current implementation offers three interaction potentials:

$$
\begin{aligned}
U(r_{i,j}) &= -\exp\Big(-\frac{r_{i,j}^2}{2\sigma^2}\Big), && \quad \text{Gauss potential, $\sigma^2=0.5$.} \\
U(r_{i,j}) &= \Big(\exp(-2ar_{i,j}) -2\exp(-ar_{i,j})\Big), && \quad \text{Morse potential, $a=2$.}\\
U(r_{i,j}) &= -\exp\Big(-\Big(\frac{r_{i,j}}{\sqrt{2\sigma^2}}\Big)^{\alpha}\Big), && \quad \text{GEM-$\alpha$ potential, $\sigma^2=0.5$, $\alpha=4$.}
\end{aligned}
$$

The system is then simulated by integrating Langevin dynamics (using one of two high-quality integrators, BAOAB [1] or UBU [2]), and observables are measured regularly. Their time series are printed to a <code>csv</code> file at the end. Optionally, the whole trajectory can be stored and printed to a file as well.  
The force calculation requires computing all pairwise distances between the particles, which is a $\mathcal{O}(N^2)$ operation---  
the computational bottleneck of the code. To speed this up, it uses multithreading via OpenMP in the force routine. At the moment, the code does
not use advanced algorithmic approaches (such as cell or Verlet lists) to improve the scaling with $N,$ since these techniques are only useful for particle counts beyond 10,000 and for interaction ranges much smaller than the simulation box. On a Dell Latitude 5530 with an i7-1265U processor, running a system of 2,000 particles for 20,000 iterations with 10 threads leads to a runtime of 1.5 minutes.

It is possible to extend the codebase by adding observables to collect, pairwise interaction functions, and Langevin dynamics integrators.
To add or change the observables to collect, modify the <code>src/measurements.h</code> file. To add or change the pairwise interaction functions, modify the 
<code>src/model.h</code> file. To add new integrators, modify the <code>src/simulation.h</code> file. The files contain comments and instructions to explain more details.

### Build From Source
To create an executable on your platform, download the files in the <code>src</code> folder and compile the `main.cpp` file. Due to the use of OpenMP, a corresponding flag will need to be passed, depending on the compiler. It might also be necessary to specify the C++17 standard (or higher). On Linux, using gcc, compilation is invoked via  
`g++ -fopenmp -O3 -std=c++17 -o simips main.cpp`.  
This will create a `simips` executable.  
If you see an error that the file `cxxopts.hpp` has not been found, download it from  
[https://github.com/jarro2783/cxxopts/tree/master/include](https://github.com/jarro2783/cxxopts/tree/master/include) and store it in your system's include folder (on Linux, this is typically /usr/include). 

### Quickstart   
To run a simulation, run the executable with `./simips`. This will simulate a single trajectory of an IPS with default parameters. It will create a `results.csv` output file holding time series data of the mean distance to the centre of mass (COM), the mean squared displacement (MSD), and the (instantaneous) kinetic temperature (Tkin) (see [3] for their definition). The first column gives the corresponding simulation times. 
The properties of the simulation can be controlled via various flags (e.g., stepsize, number of particles, number of iterations, which integrator to use etc.). Run `./simips --help` to view the possible options.

The two Python scripts to visualise the results need to be run from a Python environment with installed `numpy` and `matplotlib` packages.
To plot the time series results, run `python plot_data.py <file>`. It accepts a flag `--title plot_title` to specify a title for the plot. Run `python plot_data.py --help` for more info.

In order to create an animation of a trajectory, `./simips` needs to be executed with the `--trajectory` flag, which will prompt the program to print out trajectory data to a second file.  
This file can be read with the second Python script via `python plot_trajectory.py <file>`. Once again, there are options available (view them via the `--help` flag).

### References
[1] Leimkuhler, B., & Matthews, C. (2013). Rational construction of stochastic numerical methods for molecular sampling. *Applied Mathematics Research Express*, 2013(1), 34-56. [DOI: 10.1093/amrx/abs010](https://doi.org/10.1093/amrx/abs010) <br>
[2] Zapatero, A. Á. (2021). *Word Series for the Numerical Integration of Stochastic Differential Equations* (PhD thesis). Universidad de Valladolid. <br>
[3] Leimkuhler, B., Lohmann, R., Pavliotis, G. A., & Whalley, P. A. (2026). Cluster formation for weakly interacting kinetic Langevin dynamics. *arXiv*. [arXiv:2510.25034](https://arxiv.org/abs/2510.25034)

