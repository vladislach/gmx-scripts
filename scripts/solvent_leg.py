import os
import subprocess
import shutil

def run_command(command):
    """Helper function to run a shell command using subprocess."""
    subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def prepare_box():
    """Prepares the solvated box by running GROMACS commands."""
    # Prepare the simulation box and solvate
    run_command("gmx editconf -f output/data/merged.gro -o output/data/merged_box.gro -c -d 1.5 -bt dodecahedron")
    run_command("gmx solvate -cp output/data/merged_box.gro -cs spc216.gro -p output/data/merged.top -o output/data/merged_solvated.gro")
    run_command("mv output/data/merged.top output/data/merged_solvated.top")


def run_em():
    """Performs energy minimization for the solvated system."""
    os.makedirs("output/water/em", exist_ok=True)
    
    # Link the necessary files
    run_command("ln output/data/merged_solvated.gro output/water/em/conf.gro")
    run_command("ln output/data/merged_solvated.top output/water/em/topol.top")
    run_command("ln mdp/em.mdp output/water/em/grompp.mdp")

    # Run GROMACS energy minimization
    run_command("cd output/water/em; gmx grompp")
    run_command("cd output/water/em; gmx mdrun")


def run_nvt():
    """Performs NVT equilibration for the system."""
    os.makedirs("output/water/nvt", exist_ok=True)
    
    # Link the necessary files
    run_command("ln output/water/em/confout.gro output/water/nvt/conf.gro")
    run_command("ln output/water/em/topol.top output/water/nvt/topol.top")
    run_command("ln mdp/nvt.mdp output/water/nvt/grompp.mdp")
    
    # Run GROMACS NVT equilibration
    run_command("cd output/water/nvt; gmx grompp")
    run_command("cd output/water/nvt; gmx mdrun")


def run_npt():
    """Performs NPT equilibration for the system."""
    os.makedirs("output/water/npt", exist_ok=True)
    
    # Link the necessary files
    run_command("ln output/water/nvt/confout.gro output/water/npt/conf.gro")
    run_command("ln output/water/nvt/topol.top output/water/npt/topol.top")
    run_command("ln mdp/npt.mdp output/water/npt/grompp.mdp")
    
    # Run GROMACS NPT equilibration
    run_command("cd output/water/npt; gmx grompp")
    run_command("cd output/water/npt; gmx mdrun")


def run_fep(number_of_lambdas=11):
    """Sets up and runs the FEP calculations for the merged ligand."""
    os.makedirs("output/water/fep", exist_ok=True)

    for lamda in range(number_of_lambdas):
        lamda_dir = f"output/water/fep/lambda_{lamda}"
        os.makedirs(lamda_dir, exist_ok=True)
        
        # Link the necessary files
        run_command(f"ln output/water/npt/confout.gro {lamda_dir}/conf.gro")
        run_command(f"ln output/water/npt/topol.top {lamda_dir}/topol.top")
        
        # Modify and write the MDP file for each lambda window
        with open('mdp/fep.mdp', 'r') as f:
            mdp_content = f.read().replace("LAMBDA", str(lamda))
        with open(f"{lamda_dir}/grompp.mdp", 'w') as f:
            f.write(mdp_content)
        
        # Run GROMACS for each lambda window
        run_command(f"cd {lamda_dir}; gmx grompp -maxwarn 1")
        run_command(f"cd {lamda_dir}; gmx mdrun -dhdl dhdl.{lamda}.xvg")


def move_fep_results(number_of_lambdas=11):
    """Moves the FEP results to the analysis directory."""
    os.makedirs("output/water/analysis", exist_ok=True)

    for lamda in range(number_of_lambdas):
        src = f"output/water/fep/lambda_{lamda}/dhdl.{lamda}.xvg"
        dst = f"output/water/analysis/dhdl.{lamda}.xvg"
        if os.path.exists(src):
            shutil.move(src, dst)


# Example of running the workflow
if __name__ == "__main__":
    prepare_box()
    run_em()
    run_nvt()
    run_npt()
    run_fep()
    move_fep_results()
