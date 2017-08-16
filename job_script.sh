#!/bin/bash
#SBATCH --job-name=test_factors
#SBATCH --output=test_factors.out
#SBATCH --error=test_factors.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00

srun -n 1 R CMD BATCH ./collect_factors.R --ntasks-per-node=1
