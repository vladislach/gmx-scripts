#!/bin/bash

set -e

run_command() {
    echo "Running: $1"
    eval "$1"
}

# Step 1: Prepare Solvated Box for the Complex
prepare_box() {
    echo "Preparing solvated box for the complex..."
    mkdir -p output/complex/ions
    run_command "gmx editconf -f output/data/complex.gro -o output/data/complex_box.gro -c -d 1.5 -bt dodecahedron"
    run_command "gmx solvate -cp output/data/complex_box.gro -cs spc216.gro -p output/data/complex.top -o output/data/complex_solvated.gro"
    run_command "mv output/data/complex.top output/data/complex_solvated.top"
}

# Step 2: Add Ions to the Complex System
add_ions() {
    echo "Adding ions to neutralize the complex system..."
    run_command "ln -sf output/data/complex_solvated.gro output/complex/ions/conf.gro"
    run_command "ln -sf output/data/complex_solvated.top output/complex/ions/topol.top"
    run_command "ln -sf mdp/ions.mdp output/complex/ions/grompp.mdp"
    cd output/complex/ions
    run_command "gmx grompp -maxwarn 1"
    run_command "echo '15' | gmx genion -s topol.tpr -p -pname NA -nname CL -neutral"
    cd -
}

# Step 3: Run Energy Minimization for the Complex System
run_em() {
    echo "Running energy minimization for the complex..."
    mkdir -p output/complex/em
    run_command "ln -sf output/complex/ions/out.gro output/complex/em/conf.gro"
    run_command "ln -sf output/complex/ions/topol.top output/complex/em/topol.top"
    run_command "ln -sf mdp/em.mdp output/complex/em/grompp.mdp"
    cd output/complex/em
    run_command "gmx grompp"
    run_command "gmx mdrun"
    cd -
}

# Step 4: Run NVT Equilibration for the Complex
run_nvt() {
    echo "Running NVT equilibration for the complex..."
    mkdir -p output/complex/nvt
    run_command "ln -sf output/complex/em/confout.gro output/complex/nvt/conf.gro"
    run_command "ln -sf output/complex/em/topol.top output/complex/nvt/topol.top"
    run_command "ln -sf mdp/nvt.mdp output/complex/nvt/grompp.mdp"
    cd output/complex/nvt
    run_command "gmx grompp"
    run_command "gmx mdrun"
    cd -
}

# Step 5: Run NPT Equilibration for the Complex
run_npt() {
    echo "Running NPT equilibration for the complex..."
    mkdir -p output/complex/npt
    run_command "ln -sf output/complex/nvt/confout.gro output/complex/npt/conf.gro"
    run_command "ln -sf output/complex/nvt/topol.top output/complex/npt/topol.top"
    run_command "ln -sf mdp/npt.mdp output/complex/npt/grompp.mdp"
    cd output/complex/npt
    run_command "gmx grompp"
    run_command "gmx mdrun"
    cd -
}

# Step 6: Perform FEP Calculations for the Complex
run_fep() {
    echo "Running FEP calculations for the complex..."
    local number_of_lambdas=11
    mkdir -p output/complex/fep

    for lamda in $(seq 0 $((number_of_lambdas-1))); do
        lamda_dir="output/complex/fep/lambda_${lamda}"
        mkdir -p ${lamda_dir}
        run_command "ln -sf output/complex/npt/confout.gro ${lamda_dir}/conf.gro"
        run_command "ln -sf output/complex/npt/topol.top ${lamda_dir}/topol.top"
        
        # Modify the MDP file for each lambda window
        sed "s/LAMBDA/${lamda}/g" mdp/fep.mdp > ${lamda_dir}/grompp.mdp
        
        # Run GROMACS for each lambda window
        cd ${lamda_dir}
        run_command "gmx grompp -maxwarn 1"
        run_command "gmx mdrun -dhdl dhdl.${lamda}.xvg"
        cd -
    done
}

# Step 7: Move FEP Results to Analysis Directory
move_fep_results() {
    echo "Moving FEP results..."
    local number_of_lambdas=11
    mkdir -p output/complex/analysis

    for lamda in $(seq 0 $((number_of_lambdas-1))); do
        src="output/complex/fep/lambda_${lamda}/dhdl.${lamda}.xvg"
        dst="output/complex/analysis/dhdl.${lamda}.xvg"
        if [ -f "$src" ]; then
            mv "$src" "$dst"
        fi
    done
}

main() {
    prepare_box
    add_ions
    run_em
    run_nvt
    run_npt
    run_fep
    move_fep_results
}

main
