#!/bin/bash
#SBATCH --job-name=test_factors
#SBATCH --output=test_factors.out
#SBATCH --error=test_factors.err
#SBATCH --nodes=1
#SBATCH  --ntasks=32
#SBATCH --time=01:30:00

module load r/3.3.2-cdl


srun -n 1 --hint=nomultithread R CMD BATCH ./collect_factors.R
