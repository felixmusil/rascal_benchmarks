#!/bin/bash
#SBATCH --chdir /scratch/musil/rascal_benchmarks/scripts/spherical_invariants
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --mem 188000
#SBATCH --time 14:00:00


conda activate rascal_benchmark

python project.py run -f '{"name": "methane_sulfonic", "representation.radial_basis": "DVR"}'  --parallel 5 --progress --order random  2>&1 | tee status_41.txt
