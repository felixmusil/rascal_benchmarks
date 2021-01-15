# rascal_benchmarks

Set of benchmarks of the rascal library relying on the signac library to run and analyse the results.

## Installation

Here is the recipe to run the benchmarks from scratch:

+ create a conda environment with the required packages running
```
conda env create -f environment.yml
```
and to make the environement file run
```
conda env export > environment.yml
```

+ compile rascal (see https://github.com/cosmo-epfl/librascal for requirements) and the benchmark suite by running
```
mkdir build
cd build
cmake .. && make -j4
```

## Benchmarks

The benchmarks are composed of dedicated c++ programs (sources are in `src`) and python scripts (in `scripts`) which drive the c++ programs and also perform some benchmarks using the python interface of rascal.
Most of the parameters of the benchmarks are defined in `kernel_params.py` and `representation_params.py` and each benchmark suites selects the appropriate subsets of parameters from these file (except for the `model` and the `neighbor_list` suite).

The whole suite of benchmarks can be run with
```
cd ../scripts/
python driver.py --init --run -np 4 2>&1 | tee status.txt
```



## Datasets

The datasets have been extracted from various published articles (see `structures/raw_data/README.md` for more details) and are stored in a custom `.json` format compatible with `ASE` and `rascal` reader. These files have been produced with the `convert_structure` notebook.




rsync -avzh cosmo3:/home/musil/git/rascal_benchmarks/scripts ./

rsync -avzh --progress --max-size='5M' helvetios:/scratch/musil/rascal_benchmarks/scripts/* ./

python driver.py --run -np 15 2>&1 | tee status.txt

Sinteract -p build -c 4 -m 16G