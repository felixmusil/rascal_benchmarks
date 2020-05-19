# rascal_benchmarks


```
conda create -n rascal_benchmark python=3.7 -y
conda activate rascal_benchmark
pip install ase scikit-build
conda install -c conda-forge numpy scipy matplotlib seaborn jupyter ipython scikit-learn tqdm numba tbb mkl-service jupyter_contrib_nbextensions -y
python -m ipykernel install --user --name rascal_benchmark --display-name "rascal_benchmark"
```