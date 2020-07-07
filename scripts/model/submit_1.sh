#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/model
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --mem 180000
#SBATCH --time 48:00:00


conda activate rascal_benchmark

python project.py run -f name molecular_crystals --parallel 3 --progress  2>&1 | tee status_1.txt
