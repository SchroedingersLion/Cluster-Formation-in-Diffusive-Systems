
# Cluster Formation for weakly interacting kinetic Langevin dynamics

This folder contains information and material that allows for the reproduction of the numerical experiments for the following publication:

> **[Cluster formation for weakly interacting kinetic Langevin dynamics](https://arxiv.org/abs/2510.25034)**  
> B. Leimkuhler, R. Lohmann, G. A. Pavliotis, P. A. Whalley (in press). *Proceedings of the Royal Society A*.

### Reproducibility
The data and figures for the paper were generated using **SimIPS v1.0.0**. 
To guarantee reproducibility, do not use the `main` branch. Use the specific archived snapshot:
* **Zenodo Archive:** [https://doi.org/10.5281/zenodo.21743624](https://doi.org/10.5281/zenodo.21743624)
* **GitHub Release:** `v1.0.0 - Snapshot for publication "Cluster formation in weakly interacting kinetic Langevin dynamics"` (Git tag: `v1.0.0`)

- The particle system simulations for Sec. 4 and Appendix B in the article were run with SimIPS.

- The Jupyter notebook `stability_and_fluctuation_analysis.ipynb` implements the numerical spectral analysis of the matrix $A_{\gamma}(k,\beta)$ from the paper, as well as the generation of the corresponding plots.

- The `Animations` folder holds video animations of the particle systems (referenced in the paper).

- The animations in the `Animations` folder and the observable plots in the paper (e.g., Fig. 12) were generated with the Python scripts in `../scripts`. 


