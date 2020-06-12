# init.py
import signac
import sys
import numpy as np
from os.path import join,dirname
print(dirname(__file__))
sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import REPRESENTATION_PARAMS

project = signac.init_project('spherical_harmonics')

fixed_params = [
    'max_radial','normalize','radial_basis',
    'optimization_args','expansion_by_species_method',
    'gaussian_sigma_constant','coefficient_subselection',
    'cutoff_smooth_width','gaussian_sigma_type','soap_type','global_species',
    'inversion_symmetry','cutoff_function_parameters','cutoff_function_type',
]
fixed_params = {k:REPRESENTATION_PARAMS[0]['representation'][k] for k in fixed_params}

for REPRESENTATION_PARAM in REPRESENTATION_PARAMS:
    if np.all([REPRESENTATION_PARAM['representation'][k] == v for k,v in fixed_params.items()]):
        rep_params = {k:v for k,v in REPRESENTATION_PARAM.items() if k not in fixed_params}
        job = project.open_job(rep_params)
        job.init()