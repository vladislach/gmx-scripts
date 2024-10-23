# Relative Binding Free Energy Calculations with GROMACS

This repository contains scripts and an interactive notebook to set up and run Relative Binding Free Energy (RBFE) calculations for protein-ligand systems using **GROMACS**. The main walkthrough is provided in `example.ipynb` that demonstrates the entire process from preparing the system to performing free energy calculations with multiple lambda windows.

## Requirements

Before running any of the scripts or the notebook, ensure that GROMACS is installed, as it is required to perform molecular dynamics simulations. You can follow the official GROMACS [installation guide](https://manual.gromacs.org/current/install-guide/index.html) for instructions. It's recommended to use Conda to manage dependencies and easily set up the environment.

## Environment Setup

To set up the environment:

1. Clone the repository:

    ```bash
    git clone https://github.com/vladislach/gmx-scripts.git
    cd gmx-scripts
    ```

2. Create and activate the Conda environment:

    ```bash
    conda env create -f environment.yml
    conda activate gmx_scripts
    ```

## Notes

The main guide and walkthrough for setting up and running the RBFE calculations is provided in `example.ipynb`. Scripts for various stages of the workflow (e.g., solvation, energy minimization, equilibration, and FEP calculations) are available in the `scripts` directory.

This workflow uses [acpype](https://github.com/alanwilter/acpype) to generate ligand topologies for GROMACS and [alchemlyb](https://github.com/alchemistry/alchemlyb) to analyze the results of free energy calculations using the multistate Bennett acceptance ratio (MBAR) method. The chosen force fields are **AMBER99SB-ILDN** for the protein and **GAFF2** for the merged ligand. However, you can use any other force field compatible with GROMACS.

The `fep.mdp` file contains a sample lambda schedule for the free energy perturbation (FEP) calculations, but you can modify this file to experiment with different lambda schedules. The following resources might provide some guidance on on selecting appropriate lambda schedules:

- [Absolute Binding Free Energy from Alchemistry.org](https://alchemistry.org/wiki/Absolute_Binding_Free_Energy_-_Gromacs_2016)
- [Sample from David Mobley’s `pymbar` package](https://github.com/MobleyLab/alchemical-analysis/blob/master/samples/gromacs/inputfiles/3-methylindole-38steps/trpo.mdp)
- [Free Energy of Solvation Tutorial](https://tutorials.gromacs.org/docs/free-energy-of-solvation.html)

Additionally, to improve sampling and convergence in FEP calculations, consider using Hamiltonian replica exchange (HREX). You can enable HREX by adding the `-replex` and `-nex` flags to `mdrun` during production simulations after completing all minimization and equilibration steps. The `-multidir` flag can be used to specify the directories containing the `.tpr` files for all windows to run them in parallel. More information about these flags is available in the [GROMACS manual](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).
