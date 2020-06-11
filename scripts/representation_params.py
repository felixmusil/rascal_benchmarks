from os.path import join, abspath
from itertools import product

STRUCTURE_PATH = '../../structures'
RASCAL_BUILD_PATH = '../../build/_deps/rascal-build'
BUILD_PATH = '../../build'

# Define parameter space
fns = ['qm9.json']

misc_entries = {
    'qm9.json' : {
        "N_ITERATIONS" : 100,
        "n_structures" : 50,
    }
}

# spherical expansion params
ns = [3,4,5]
ls = [3]
grads = [False]
rcs = [3]
gaussian_sigmas = [0.3]
normalize = [False]
soap_types = ['PowerSpectrum']
radial_basis = ['GTO']
expansion_by_species_method = ['structure wise']
feature_sparsification_fractions = [None]

optimization_args = {
    'GTO' : {},
}

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
                 radial_basis=rb, normalize=norm,
                 optimization_args=optimization_args[rb],
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
    sp.update(misc_entries[fn])
    REPRESENTATION_PARAMS.append(sp)