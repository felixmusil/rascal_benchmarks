#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_invariants
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --mem 188000
#SBATCH --time 06:00:00


conda activate rascal_benchmark

python project.py run -f name methane_liquid --parallel 5  --order random  2>&1 | tee status_3.txt
