#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/radial_integral
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 30
#SBATCH --mem 180000
#SBATCH --time 10:00:00


conda activate rascal_benchmark

python project.py run --parallel 30
