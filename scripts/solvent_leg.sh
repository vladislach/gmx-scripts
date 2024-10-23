#!/bin/bash

set -e

run_command() {
    echo "Running: $1"
    eval "$1"
}

# Step 1: Prepare Solvated Box
prepare_box() {
    echo "Preparing solvated box..."
    run_command "gmx editconf -f output/data/merged.gro -o output/data/merged_box.gro -c -d 1.5 -bt dodecahedron"
    run_command "gmx solvate -cp output/data/merged_box.gro -cs spc216.gro -p output/data/merged.top -o output/data/merged_solvated.gro"
    run_command "mv output/data/merged.top output/data/merged_solvated.top"
}

# Step 2: Run Energy Minimization
run_em() {
    echo "Running energy minimization..."
    mkdir -p output/water/em
    run_command "ln -sf output/data/merged_solvated.gro output/water/em/conf.gro"
    run_command "ln -sf output/data/merged_solvated.top output/water/em/topol.top"
    run_command "ln -sf mdp/em.mdp output/water/em/grompp.mdp"
    cd output/water/em
    run_command "gmx grompp"
    run_command "gmx mdrun"
    cd -
}

# Step 3: Run NVT Equilibration
run_nvt() {
    echo "Running NVT equilibration..."
    mkdir -p output/water/nvt
    run_command "ln -sf output/water/em/confout.gro output/water/nvt/conf.gro"
    run_command "ln -sf output/water/em/topol.top output/water/nvt/topol.top"
    run_command "ln -sf mdp/nvt.mdp output/water/nvt/grompp.mdp"
    cd output/water/nvt
    run_command "gmx grompp"
    run_command "gmx mdrun"
    cd -
}

# Step 4: Run NPT Equilibration
run_npt() {
    echo "Running NPT equilibration..."
    mkdir -p output/water/npt
    run_command "ln -sf output/water/nvt/confout.gro output/water/npt/conf.gro"
    run_command "ln -sf output/water/nvt/topol.top output/water/npt/topol.top"
    run_command "ln -sf mdp/npt.mdp output/water/npt/grompp.mdp"
    cd output/water/npt
    run_command "gmx grompp"
    run_command "gmx mdrun"
    cd -
}

# Step 5: Run FEP for Multiple Lambda Windows
run_fep() {
    echo "Running FEP calculations..."
    local number_of_lambdas=11
    mkdir -p output/water/fep

    for lamda in $(seq 0 $((number_of_lambdas-1))); do
        lamda_dir="output/water/fep/lambda_${lamda}"
        mkdir -p ${lamda_dir}
        run_command "ln -sf output/water/npt/confout.gro ${lamda_dir}/conf.gro"
        run_command "ln -sf output/water/npt/topol.top ${lamda_dir}/topol.top"
        
        # Modify the MDP file for each lambda window
        sed "s/LAMBDA/${lamda}/g" mdp/fep.mdp > ${lamda_dir}/grompp.mdp
        
        # Run GROMACS for each lambda window
        cd ${lamda_dir}
        run_command "gmx grompp -maxwarn 1"
        run_command "gmx mdrun -dhdl dhdl.${lamda}.xvg"
        cd -
    done
}

# Step 6: Move FEP Results to Analysis Directory
move_fep_results() {
    echo "Moving FEP results..."
    local number_of_lambdas=11
    mkdir -p output/water/analysis

    for lamda in $(seq 0 $((number_of_lambdas-1))); do
        src="output/water/fep/lambda_${lamda}/dhdl.${lamda}.xvg"
        dst="output/water/analysis/dhdl.${lamda}.xvg"
        if [ -f "$src" ]; then
            mv "$src" "$dst"
        fi
    done
}

main() {
    prepare_box
    run_em
    run_nvt
    run_npt
    run_fep
    move_fep_results
}

main
