#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_invariants
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --mem 188000
#SBATCH --time 12:00:00


conda activate rascal_benchmark

python project.py run -f name qm9 --parallel 3 --progress --order random  2>&1 | tee status_0.txt
