# init.py
import signac
import sys
from copy import deepcopy
import numpy as np
from os.path import join,dirname

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
del_params = [
    'normalize','radial_basis',
    'optimization_args','expansion_by_species_method',
    'gaussian_sigma_constant','coefficient_subselection',
    'soap_type','global_species',
    'inversion_symmetry','cutoff_function_parameters','cutoff_function_type',
]

fixed_params = {k:REPRESENTATION_PARAMS[0]['representation'][k] for k in fixed_params}

for REPRESENTATION_PARAM in REPRESENTATION_PARAMS:
    if REPRESENTATION_PARAM['nl_type'] != 'Full':
        pass
    if np.all([REPRESENTATION_PARAM['representation'][k] == v for k,v in fixed_params.items()]):
        rep_args = deepcopy(REPRESENTATION_PARAM)
        rep_args['N_ITERATIONS'] = 50
        rep_args['representation'] = {
            k:v for k,v in rep_args['representation'].items() if k not in del_params}
        job = project.open_job(rep_args)
        job.init()