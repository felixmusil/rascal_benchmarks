# init.py
import signac
import sys
from copy import deepcopy
import numpy as np
from os.path import join,dirname
sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH

project = signac.init_project('model')
# Define parameter space

names = ['silicon_bulk']

fns = {
    'silicon_bulk': 'silicon_bulk.json',
    'methane_liquid': 'methane_liquid.json',
    'methane_sulfonic': 'methane_sulfonic.json'
}

misc_entries = {
    'silicon_bulk' : {
        "N_ITERATIONS" : 20,
        'start_structure': 600,
        "n_structures" : 20,
    },
    'methane_liquid' : {
        "N_ITERATIONS" : 5,
        'start_structure': 100,
        "n_structures" : 35,
    },
    'methane_sulfonic' : {
        "N_ITERATIONS" : 5,
        'start_structure': 0,
        "n_structures" : 85,
    },
}

global_species = {
    'qm9' : [1,6,7,8,9],
    'molecular_crystals' : [1, 6, 7, 8],
    'silicon_bulk' : [14],
    'methane_liquid' : [1, 6],
    'methane_sulfonic' : [1,6,8,16],
}
# soap l_max=12 n_max=10 atom_sigma=0.5 zeta=4 cutoff=5.0 cutoff_transition_width=1.0 central_weight=1.0 n_sparse=9000 delta=3.0 f0=0.0 covariance_type=dot_product sparse_method=cur_points
# signal_variance="3.0000000000000000" signal_mean=".00000000000000000E+000" sparsified="T" n_permutations="1" covariance_type="2" zeta="4.
models = {
    'silicon_bulk' : [
        {
         'representation' :
            dict(
            interaction_cutoff=5., cutoff_smooth_width=1.,
            max_radial=10, max_angular=12, gaussian_sigma_type="Constant",
            soap_type="PowerSpectrum",
            normalize=True,
            expansion_by_species_method='structure wise',
            global_species=global_species['silicon_bulk'],
            compute_gradients=False,
            cutoff_function_parameters=dict(),
            cutoff_function_type="ShiftedCosine",
            gaussian_sigma_constant=0.5,
            coefficient_subselection=None,
            ),
         'kernel':
           dict(
               name='GAP', zeta=4, target_type='Structure', kernel_type='Sparse'
           ),
         'RandomFilter': dict(
             Nselect={14:100}, act_on='sample per species'
         ),

        },
    ],
}

self_contributions = {
    'silicon_bulk' : {14: -158.54496821},
}


for name in names:
    for model in models[name]:
        rep_args = deepcopy(model)
        rep_args['name'] = name
        rep_args['filename'] = join(STRUCTURE_PATH,fns[name])
        rep_args['self_contributions'] = self_contributions[name]
        rep_args['train_with_grad'] = False
        rep_args.update(misc_entries[name])
        job = project.open_job(rep_args)
        job.init()