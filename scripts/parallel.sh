#!/bin/bash
#SBATCH --job-name=gmx_mdrun_multidir
#SBATCH --output=gmx_mdrun_multidir_%j.log
#SBATCH --error=gmx_mdrun_multidir_%j.err
#SBATCH --ntasks=16                      # Number of MPI tasks (change as needed)
#SBATCH --cpus-per-task=4                # Number of OpenMP threads per task (change as needed)
#SBATCH --time=48:00:00                  # Time limit hrs:min:sec
#SBATCH --partition=gpu                  # Specify the partition (queue) to submit to
#SBATCH --mem=8G                         # Memory per node (adjust based on requirements)

# Load the GROMACS module (depends on your cluster environment)
module load GROMACS

# Run mdrun with the -multidir option, distributing tasks across the directories
# The -multidir flag enables running simulations in multiple directories simultaneously.
srun gmx mdrun -s topol.tpr -deffnm prod -replex 1000 \
    -multidir output/water/fep/lambda_* -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK