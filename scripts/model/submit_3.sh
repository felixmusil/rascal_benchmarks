#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 3
#SBATCH --mem 188000
#SBATCH --time 24:00:00


conda activate rascal_benchmark
cd /scratch/musil/rascal_benchmarks/scripts/model
python project.py run -f name methane_liquid --parallel 3 --progress --order random  2>&1 | tee status_3.txt
