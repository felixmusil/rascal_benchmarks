from os.path import join, abspath
from itertools import product
from copy import deepcopy
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH

# Define parameter space
fns = {
    'qm9': 'qm9.json',
    'molecular_crystals': 'molecular_crystals_50.json',
    'silicon_bulk': 'silicon_bulk.json',
    'methane_liquid': 'methane_liquid.json',
    'methane_sulfonic': 'methane_sulfonic_150.json'
}

misc_entries = {
    'qm9' : {
        "N_ITERATIONS" : 20,
        "n_structures" : 500,
        'start_structure': 0,
    },
    'molecular_crystals' : {
        "N_ITERATIONS" : 20,
        'start_structure': 0,
        "n_structures" : 30,
    },
    'silicon_bulk' : {
        "N_ITERATIONS" : 20,
        'start_structure': 600,
        "n_structures" : 100,
    },
    'methane_liquid' : {
        "N_ITERATIONS" : 20,
        'start_structure': 100,
        "n_structures" : 50,
    },
    'methane_sulfonic' : {
        "N_ITERATIONS" : 20,
        'start_structure': 0,
        "n_structures" : 150,
    },
}

# spherical expansion params
ns = range(2,16,2)
ls = range(1,16,2)
grads = [False, True]
rcs = [4]
gaussian_sigmas = [0.3]
normalize = [False, True]
soap_types = ['PowerSpectrum']
radial_basis = ['GTO','DVR','GTO_Spline','DVR_Spline']
expansion_by_species_method = ['structure wise']
feature_sparsification_fractions = [None,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6]
seed = 10

radial_basis_args = {
    'GTO' : {'radial_basis': 'GTO', 'optimization_args': {}},
    'DVR' : {'radial_basis': 'DVR', 'optimization_args': {}},
    'GTO_Spline' : {
            'radial_basis': 'GTO',
            'optimization_args':{
                'type': 'Spline', 'accuracy': 1e-8, 'range': [0,0]
            }
    },
    'DVR_Spline' : {
            'radial_basis': 'DVR',
            'optimization_args':{
                'type': 'Spline', 'accuracy': 1e-8, 'range': [1e-6,0]
            }
    }
}

accuracies = [1e-2,1e-4,1e-6,1e-8,1e-10]

global_species = {
    'qm9' : [1,6,7,8,9],
    'molecular_crystals' : [1, 6, 7, 8],
    'silicon_bulk' : [14],
    'methane_liquid' : [1, 6],
    'methane_sulfonic' : [1,6,8,16],
}

neighbourlist = {
    "Full": [
        {"initialization_arguments": {"cutoff": -1}, "name": "neighbourlist"},
        {"initialization_arguments": {}, "name": "centercontribution"},
        {"initialization_arguments": {"cutoff": -1}, "name": "strict"}
    ],
    "Half": [
        {"initialization_arguments": {"cutoff": -1}, "name": "neighbourlist"},
        {"initialization_arguments": {}, "name": "halflist"},
        {"initialization_arguments": {}, "name": "centercontribution"},
        {"initialization_arguments": {"cutoff": -1}, "name": "strict"}
    ],
}


REPRESENTATION_PARAMS = []

for ((name, fn), n, l, grad, rc, gaussian_sigma, norm, soap_type, rb,
    exp_meth, feature_sparsification_fraction, (nl_type, nl)) in product(
        fns.items(), ns, ls, grads, rcs, gaussian_sigmas, normalize, soap_types, radial_basis,
        expansion_by_species_method, feature_sparsification_fractions,
        neighbourlist.items()):
    adaptors = deepcopy(nl)
    for el in adaptors:
        if 'cutoff' in el["initialization_arguments"]:
            el["initialization_arguments"]['cutoff'] = rc
    sp = {
        'name' : name,
        'filename': join(STRUCTURE_PATH,fn),
        'representation': dict(
            interaction_cutoff=rc, cutoff_smooth_width=0.5,
            max_radial=n, max_angular=l, gaussian_sigma_type="Constant",
            soap_type=soap_type,
            normalize=norm,
            expansion_by_species_method=exp_meth,
            global_species=global_species[name],
            compute_gradients=grad,
            inversion_symmetry=True,
            cutoff_function_parameters=dict(),
            cutoff_function_type="ShiftedCosine",
            gaussian_sigma_constant=gaussian_sigma,
            coefficient_subselection=feature_sparsification_fraction,
        ),
        'nl_type' : nl_type,
        'adaptors' : adaptors,
        'seed':seed,
    }
    if rb not in ['GTO', 'DVR']:
        kargs = deepcopy(radial_basis_args[rb])
        for accuracy in accuracies:
            sp_ = deepcopy(sp)
            kargs_ = deepcopy(kargs)
            kargs_['optimization_args']['accuracy'] = accuracy
            kargs_['optimization_args']['range'][1] = rc
            sp_['representation'].update(**kargs_)
            sp_.update(misc_entries[name])
            REPRESENTATION_PARAMS.append(sp_)
    else:
        kargs = deepcopy(radial_basis_args[rb])
        sp['representation'].update(**kargs)
        sp.update(misc_entries[name])
        REPRESENTATION_PARAMS.append(sp)