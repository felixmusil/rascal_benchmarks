# init.py
import signac
import sys
import numpy as np
from copy import deepcopy
from os.path import join,dirname

sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import REPRESENTATION_PARAMS

project = signac.init_project('radial_integral')

fixed_params = [
    'normalize','expansion_by_species_method',
    'coefficient_subselection','global_species',
    'cutoff_smooth_width','soap_type',
    'inversion_symmetry','cutoff_function_parameters','cutoff_function_type',
]
del_params = [
    'normalize','expansion_by_species_method',
    'coefficient_subselection','global_species',
    'soap_type',
    'inversion_symmetry','cutoff_function_parameters','cutoff_function_type',
]

fixed_params = {k:REPRESENTATION_PARAMS[0]['representation'][k] for k in fixed_params}

for REPRESENTATION_PARAM in REPRESENTATION_PARAMS:
    if REPRESENTATION_PARAM['nl_type'] != 'Full':
        continue
    if REPRESENTATION_PARAM['name'] != 'qm9':
        continue

    if np.all([REPRESENTATION_PARAM['representation'][k] == v for k,v in fixed_params.items()]):
        rep_args = deepcopy(REPRESENTATION_PARAM)
        rep_args['N_ITERATIONS'] = 50
        rep_args['n_structures'] = 2500
        rep_args['representation'] = {
            k:v for k,v in rep_args['representation'].items() if k not in del_params}
        job = project.open_job(rep_args)
        job.init()