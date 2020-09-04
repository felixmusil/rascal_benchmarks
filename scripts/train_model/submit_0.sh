#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 188000
#SBATCH --time 12:00:00

cd /scratch/musil/rascal_benchmarks/scripts/train_model

conda activate rascal_benchmark

python project.py run --parallel 2 --progress --order random  2>&1 | tee status_0.txt
