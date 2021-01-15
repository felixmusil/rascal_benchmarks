# init.py
import signac
import sys
from copy import deepcopy
import numpy as np
from os.path import join,dirname
from itertools import product
sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH

project = signac.init_project('train_model')
# Define parameter space

names = ['silicon_bulk','methane_liquid']
names = ['methane_liquid']
seed = 10
fns = {
    'silicon_bulk': 'silicon_bulk.json',
    'methane_liquid': 'methane_liquid.json',
    'molecular_crystals': 'molecular_crystals_100.json',
    'qm9': 'qm9.json',
    'methane_sulfonic': 'methane_sulfonic_150.json',
}

misc_entries = {
    'silicon_bulk' : {
        'start_structure': 1,
         "n_structures" : 2474,
        # "n_structures" : 200,
    },
    'methane_liquid' : {
        'start_structure': 1,
         "n_structures" : 20,
        # "n_structures" : 258,
    },
}


global_species = {
    'qm9' : [1,6,7,8,9],
    'molecular_crystals' : [1, 6, 7, 8],
    'silicon_bulk' : [14],
    'methane_liquid' : [1, 6],
    'methane_sulfonic' : [1,6,8,16],
}

n_sparse = [100, 200, 500, 1000, 2000, 5000, 7000, 9000]
# n_sparse = [200]
sparse_point_subselections = {
    'silicon_bulk' : [{14:v} for v in n_sparse],
    'methane_liquid' : [{1:int(v*0.8),6:int(v*0.2)} for v in n_sparse],
}

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
            radial_basis="GTO",
            optimization_args={
                  "type": "Spline", "accuracy": 1e-08, "range": [0, 5]
                },
            ),
         'kernel':
           dict(
               name='GAP', zeta=4, target_type='Structure', kernel_type='Sparse'
           ),
        },
    ],
    'methane_liquid' : [
        {
         'representation' :
            dict(
            interaction_cutoff=5, cutoff_smooth_width=1.,
            max_radial=8, max_angular=6, gaussian_sigma_type="Constant",
            soap_type="PowerSpectrum",
            normalize=True,
            expansion_by_species_method='structure wise',
            global_species=global_species['methane_liquid'],
            compute_gradients=False,
            cutoff_function_parameters=dict(),
            cutoff_function_type="ShiftedCosine",
            gaussian_sigma_constant=0.4,
            coefficient_subselection=None,
            radial_basis="GTO",
            optimization_args={
                  "type": "Spline", "accuracy": 1e-08, "range": [0, 5]
                },
            ),
         'kernel':
           dict(
               name='GAP', zeta=2, target_type='Structure', kernel_type='Sparse'
           ),
        },
    ],
}

f_feature = [0.02, 0.05, 0.1, 0.2,0.4, 0.5, 0.7, 1]

feature_subselections = {}
for name,dd in models.items():
    aa = dd[0]['representation']
    # number of PowerSpectrum features in QUIP
    n_feat = int((aa['max_angular']+1) * aa['max_radial']*len(global_species[name])*(aa['max_radial']*len(global_species[name])+1) /2)
    feature_subselections[name] = [int(v*n_feat) for v in f_feature]

self_contributions = {
    'silicon_bulk' : {14: -158.54496821},
    'qm9': {1: 0, 6: 0, 7: 0, 8: 0, 9: 0},
    'molecular_crystals': {1: -1.2066278536666357, 6: -18.42141665763611, 7: -28.21055160856482, 8: -41.63852285569995},
    'methane_liquid': {1: 0, 6: 0},
    'methane_sulfonic': {1: -0.6645519125911715, 6: -5.654232251386078, 8: -15.852522852103935, 16: -9.17258361289801}
}
if __name__ == '__main__':
    for name in names:
        for model, sparse_point_subselection, feature_subselection in product(models[name], sparse_point_subselections[name], feature_subselections[name]):
            rep_args = deepcopy(model)
            rep_args['name'] = name
            rep_args['filename'] = join(STRUCTURE_PATH,fns[name])
            rep_args['self_contributions'] = self_contributions[name]
            rep_args['train_with_grad'] = True
            rep_args['sparse_point_subselection'] = sparse_point_subselection
            rep_args['feature_subselection'] = feature_subselection
            rep_args.update(misc_entries[name])
            job = project.open_job(rep_args)
            job.init()