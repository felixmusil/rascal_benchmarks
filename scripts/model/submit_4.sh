#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 189000
#SBATCH --time 06:00:00


conda activate rascal_benchmark
cd /scratch/musil/rascal_benchmarks/scripts/model
python project.py run -f name methane_sulfonic --parallel 10 --progress --order random  2>&1 | tee status_4.txt
