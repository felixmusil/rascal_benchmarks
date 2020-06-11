from itertools import product
from .representation_params import fns, REPRESENTATION_PARAMS

n_sparses = {
    'qm9.json' : {1:10, 6:10, 7:10, 8:10},
}


target_type = {
    'qm9.json' : 'Structure'
}

kernel_types = ['Sparse']
zetas = [1]
kernel_name = {
    'Sparse' : 'GAP',
}

KERNEL_PARAMS = []
for (REPRESENTATION_PARAM, n_sparse, kernel_type, zeta) in product(REPRESENTATION_PARAMS, fns, n_sparses, kernel_types, zetas):
    fn = REPRESENTATION_PARAM['filename']
    sp = {
        'kernel' : dict(
                 name=kernel_name[kernel_type], kernel_type=kernel_type,
                 target_type=target_type[fn], zeta = zeta,
        )
    }
    KERNEL_PARAMS.append(sp)
