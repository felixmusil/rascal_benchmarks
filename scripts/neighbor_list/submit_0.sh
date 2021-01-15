#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 36
#SBATCH --mem 180000
#SBATCH --time 01:00:00
#SBATCH -p debug

conda activate rascal_benchmark
cd /scratch/musil/rascal_benchmarks/scripts/neighbor_list
python project.py run  --parallel 18 --progress --order random 2>&1 | tee status_0.txt
