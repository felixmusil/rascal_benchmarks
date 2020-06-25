#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_harmonics
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 180000
#SBATCH --time 10:00:00


conda activate rascal_benchmark

python project.py run  --parallel 20
