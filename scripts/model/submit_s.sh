#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/model
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 188000
#SBATCH --time 24:00:00


conda activate rascal_benchmark

python project.py run -j 4d01adbaf9257e624c24473bdd1826e3 113a03b32340e818fe985ec976999bfc 0c715d54941df61731e91bed7832d10f b436e92c30a834020747b2a9ff5decc0 ff9d74ea4df4edf362df5f9f1d415f62  deabe4d997075016f3e8973b3e3c5345 7e93a36dadcb68b88cf11131c0504b8c 7608d80ce2206b98fcb51eee975a2cfe 9931702bf91d6222bf9a3292071cc5bd  --parallel 3 --progress 2>&1 | tee status_s.txt
