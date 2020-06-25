#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_expansion
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 36
#SBATCH --mem 188096
#SBATCH --time 48:00:00


conda activate rascal_benchmark

python project.py run --parallel 20 2>&1 | tee status.txt
