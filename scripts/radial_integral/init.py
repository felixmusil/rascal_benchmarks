# init.py
import signac
import sys
import numpy as np
from os.path import join,dirname
print(dirname(__file__))
sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import REPRESENTATION_PARAMS

project = signac.init_project('radial_integral')

fixed_params = [
    'normalize','expansion_by_species_method',
    'coefficient_subselection','global_species',
    'cutoff_smooth_width','expansion_by_species_method',
    'inversion_symmetry','cutoff_function_parameters','cutoff_function_type',
]
fixed_params = {k:REPRESENTATION_PARAMS[0]['representation'][k] for k in fixed_params}

for REPRESENTATION_PARAM in REPRESENTATION_PARAMS:
    if np.all([REPRESENTATION_PARAM['representation'][k] == v for k,v in fixed_params.items()]):
        rep_params = {k:v for k,v in REPRESENTATION_PARAM.items() if k not in fixed_params}
        job = project.open_job(REPRESENTATION_PARAM)
        job.init()