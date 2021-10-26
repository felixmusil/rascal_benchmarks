rascal_benchmarks
-----------------

Set of benchmarks of the rascal library relying on the signac library to run and analyse the results.
These benchmarks have been used to produce the results for the journal article [Efficient implementation of atom-density representations. J. Chem. Phys. **154**, 114109 (2021)](https://doi.org/10.1063/5.0044689).

Installation
============

Here is the recipe to run the benchmarks from scratch:

+ create a conda environment with the required python packages running
```
conda env create -f environment.yml
```
and to make the environement file run
```
conda env export > environment.yml
```

+ compile rascal (see https://github.com/cosmo-epfl/librascal for c++ requirements) and the benchmark suite by running
```
mkdir build
cd build
cmake .. && make -j4
```

For the QUIP benchmarks, [QUIP](https://github.com/libatoms/QUIP) is required
and the `quip` and `gap_fit` binaries should be in the `$PATH`.

Benchmarks
==========

The benchmarks are composed of dedicated c++ programs (sources are in `src`) and python scripts (in `scripts`) which drive the c++ programs and also perform some benchmarks using the python interface of rascal.
Most of the parameters of the benchmarks are defined in `kernel_params.py` and `representation_params.py` and each benchmark suites selects the appropriate subsets of parameters from these file (except for the `model` and the `neighbor_list` suite).

The whole suite of librascal benchmarks can be run with
```
cd ../scripts/
python driver.py --init --run -np 4 2>&1 | tee status.txt
```
or on a slurm compatible HPC
```
cd ../scripts/
python driver.py --init --submit
```
where the provided `submit_*.sh` where used on the helvetios@EPFL.

The QUIP benchmarks can be run using `signac` in the `scripts/quip` directory.
The complete process of running the benchmarks involves first training potentials
for the different datasets and parameter combinations, then evaluating those
potentials to get timings.  Most of this workflow is automated using `signac`,
so a `python project.py submit` should be sufficient to run each step of the
process (use `python project.py status` to inspect execution progress), given
adequate computing resources and an appropriate
[signac cluster configuration](https://docs.signac.io/en/latest/cluster_submission.html).

Datasets
========

The datasets have been extracted from various published articles (see `structures/raw_data/README.md` for more details) and are stored in a custom `.json` format compatible with `ASE` and `rascal` reader. These files have been produced with the `convert_structure` notebook.


Notebooks
=========

Some analysis of the benchmarks results have been done with the `model_prediction_timings` and `query4plot_structure2rep` notebooks. The silicon MLIPs used to predict structural properties have been developed in the `MLIP_Si-prop` notebook.  Analysis and plotting of the QUIP benchmarks was done in the `quip_timings.ipynb` notebook.

some useful commands to work remotely
=====================================

rsync -avzh cosmo3:/home/musil/git/rascal_benchmarks/scripts ./

rsync -avzh --progress --max-size='5M' helvetios:/scratch/musil/rascal_benchmarks/scripts/* ./

python driver.py --run -np 15 2>&1 | tee status.txt

Sinteract -p build -c 4 -m 16G
