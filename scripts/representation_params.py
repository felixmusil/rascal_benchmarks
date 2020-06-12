from os.path import join, abspath
from itertools import product
from copy import deepcopy

STRUCTURE_PATH = '../../structures'
RASCAL_BUILD_PATH = '../../build/_deps/rascal-build'
BUILD_PATH = '../../build'

# Define parameter space
fns = ['qm9.json']

misc_entries = {
    'qm9.json' : {
        "N_ITERATIONS" : 20,
        "n_structures" : 500,
    }
}

# spherical expansion params
ns = range(2,16,2)
ls = range(1,16,2)
grads = [False]
rcs = [2,3,4,5]
gaussian_sigmas = [0.3]
normalize = [False, True]
soap_types = ['PowerSpectrum']
radial_basis = ['GTO','DVR','GTO_Spline','DVR_Spline']
expansion_by_species_method = ['structure wise']
feature_sparsification_fractions = [None]

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
                'type': 'Spline', 'accuracy': 1e-8, 'range': [0,0]
            }
    }
}

accuracies = [1e-2,1e-4,1e-6,1e-8,1e-10]

global_species = {
    'qm9.json' : [1,6,7,8],
}

REPRESENTATION_PARAMS = []

for (fn, n, l, grad, rc, gaussian_sigma, norm, soap_type, rb,
    exp_meth, feature_sparsification_fraction) in product(
        fns, ns, ls, grads, rcs, gaussian_sigmas, normalize, soap_types, radial_basis,
        expansion_by_species_method, feature_sparsification_fractions):

    sp = {
        'filename': join(STRUCTURE_PATH,fn),
        'representation': dict(
            interaction_cutoff=rc, cutoff_smooth_width=0.5,
                 max_radial=n, max_angular=l, gaussian_sigma_type="Constant",
                 soap_type=soap_type,
                 normalize=norm,
                 expansion_by_species_method=exp_meth,
                 global_species=global_species[fn],
                 compute_gradients=grad,
                 inversion_symmetry=True,
                 cutoff_function_parameters=dict(),
                 cutoff_function_type="ShiftedCosine",
                 gaussian_sigma_constant=gaussian_sigma,
                 coefficient_subselection=feature_sparsification_fraction
        )
    }
    if rb not in ['GTO', 'DVR']:
        kargs = deepcopy(radial_basis_args[rb])
        for accuracy in accuracies:
            kargs['optimization_args']['accuracy'] = accuracy
            kargs['optimization_args']['range'][1] = rc
            sp['representation'].update(**kargs)
            sp.update(misc_entries[fn])
            REPRESENTATION_PARAMS.append(sp)
    else:
        kargs = deepcopy(radial_basis_args[rb])
        sp['representation'].update(**kargs)
        sp.update(misc_entries[fn])
        REPRESENTATION_PARAMS.append(sp)