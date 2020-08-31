#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/model
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --mem 189000
#SBATCH --time 01:00:00


conda activate rascal_benchmark

python project.py run -f name methane_sulfonic --parallel 3 --progress --order random  2>&1 | tee status_4.txt
