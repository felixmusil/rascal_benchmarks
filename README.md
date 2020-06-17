# rascal_benchmarks



```
conda env export > environment.yml
```

```
conda env create -f environment.yml
```


rsync -avzh cosmo3:/home/musil/git/rascal_benchmarks/scripts ./

python driver.py --run -np 15 2>&1 | tee status.txt