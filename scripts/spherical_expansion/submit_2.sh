#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_expansion
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 11
#SBATCH --mem 189000
#SBATCH --time 06:00:00


conda activate rascal_benchmark

python project.py run -f name silicon_bulk --parallel 10 --progress --order random  2>&1 | tee status_2.txt
