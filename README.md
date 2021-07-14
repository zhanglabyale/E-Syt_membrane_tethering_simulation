# E-Syt membrane tethering simulation
The MATLAB codes here are used to calculate the probabilities, forces, and equilibrium energies of transmembrane binding of a single E-Syt1 or E-Syt2 molecule described in the following paper

Ge, J. et al. Stepwise membrane binding of the C2 domains of the extended synaptotagmins revealed by optical tweezers. Nat. Chem. Biol. Under peer review. https://www.researchsquare.com/article/rs-523346/v1 (2021).

Running the main code Esyts_tether_simulation.m will output a figure containing plots shown in Figs. 6b, 6d, and Supplementary Figs. 9b, & 9d of the paper, which generally takes less than 10 seconds in a desktop computer. An example output is shown in Output_figure.fig.

Other MATLAB codes (conc_tethered_fjc.m, force_energy2.m, membrane_potential_f2.m, MS.m, and pullforce1.m) are functions required by Esyts_tether_simulation.m.

The codes are tested on MATLAB version 9.4.0.813654 (R2018a).
