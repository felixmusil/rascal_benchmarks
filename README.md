# rascal_benchmarks



```
conda env export > environment.yml
```

```
conda env create -f environment.yml
```

```
mkdir build
cd build
cmake .. && make -j4
cd ../scripts/
python driver.py --init --run -np 4 2>&1 | tee status.txt
```

rsync -avzh cosmo3:/home/musil/git/rascal_benchmarks/scripts ./

rsync -avzh helvetios:/scratch/musil/rascal_benchmarks/scripts ./

python driver.py --run -np 15 2>&1 | tee status.txt
