# init.py
import signac
import sys
from copy import deepcopy
import numpy as np
from os.path import join,dirname

sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import REPRESENTATION_PARAMS

project = signac.init_project('spherical_invariants')

fixed_params = []
del_params = fixed_params
fixed_params = {k:REPRESENTATION_PARAMS[0]['representation'][k] for k in fixed_params}
fixed_accuracy = 1e-8

misc_entries = {
    'qm9' : {
        "N_ITERATIONS" : 5,
        "n_structures" : 300,
        'start_structure': 0,
    },
    'molecular_crystals' : {
        "N_ITERATIONS" : 5,
        'start_structure': 0,
        "n_structures" : 20,
    },
    'silicon_bulk' : {
        "N_ITERATIONS" : 5,
        'start_structure': 600,
        "n_structures" : 60,
    },
    'methane_liquid' : {
        "N_ITERATIONS" : 5,
        'start_structure': 100,
        "n_structures" : 35,
    },
    'methane_sulfonic' : {
        "N_ITERATIONS" : 5,
        'start_structure': 0,
        "n_structures" : 10,
    },
}


for REPRESENTATION_PARAM in REPRESENTATION_PARAMS:
    if np.all([REPRESENTATION_PARAM['representation'][k] == v for k,v in fixed_params.items()]):
        rpr = REPRESENTATION_PARAM['representation']
        if 'accuracy' in rpr['optimization_args']:
            if rpr['optimization_args']['accuracy'] != fixed_accuracy:
                continue
        # avoid sparsification when too many or to few features
        if rpr['coefficient_subselection'] is not None:
            if rpr['max_radial'] < 6 or rpr['max_angular'] < 6:
                continue
            if rpr['max_radial'] > 12 or rpr['max_angular'] > 12:
                continue
            if rpr['radial_basis'] != 'GTO' or 'accuracy' not in rpr['optimization_args']:
                continue

        rep_args = deepcopy(REPRESENTATION_PARAM)

        rep_args.update(**misc_entries[rep_args['name']])
        rep_args['representation'] = {
            k:v for k,v in rep_args['representation'].items() if k not in del_params}
        job = project.open_job(rep_args)
        job.init()