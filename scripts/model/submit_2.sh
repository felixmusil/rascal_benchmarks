#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --mem 188000
#SBATCH --time 24:00:00


conda activate rascal_benchmark
cd /scratch/musil/rascal_benchmarks/scripts/model
python project.py run -f name silicon_bulk --parallel 3 --progress --progress  --order random 2>&1 | tee status_2.txt
