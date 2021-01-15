#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 36
#SBATCH --mem 180000
#SBATCH --time 03:00:00


conda activate rascal_benchmark
cd /scratch/musil/rascal_benchmarks/scripts/radial_integral
python project.py run --parallel 26 --progress --order random 2>&1 | tee status_0.txt
