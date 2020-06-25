#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_expansion
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 180000
#SBATCH --time 48:00:00


conda activate rascal_benchmark

python project.py run -f name molecular_crystals --parallel 10
