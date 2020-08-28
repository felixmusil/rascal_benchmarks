"""Initialize the data space for QUIP overall model benchmarks"""

import itertools
import os

import signac

project = signac.init_project('quip-benchmarks')

system_names = [
    'silicon_bulk',
    'methane_liquid',
    'methane_sulfonic',
    'molecular_crystals',
    'qm9'
]

system_filenames = {
    'silicon_bulk': 'gp_iter6_sparse9k.xml.xyz',
    'methane_liquid': 'repo-fit-bulk/mebox-minimal-pbe0-b1b.xyz',
    'methane_sulfonic': 'methane_sulfonic_150.xyz',
    'molecular_crystals': 'CSD1000-r.xyz',
    'qm9': 'qm9_5000_dHf_peratom.xyz'
}

system_energy_keys = {
    'silicon_bulk': 'dft_energy',
    'methane_liquid': 'energy',
    'methane_sulfonic': 'energy',
    'molecular_crystals': 'energy',
    'qm9': 'energy'
}

system_force_keys = {
    'silicon_bulk': 'dft_force',
    'methane_liquid': 'force',
    'methane_sulfonic': 'force',
    'molecular_crystals': 'none',
    'qm9': 'none'
}

global_species = {
    'qm9' : [1,6,7,8,9],
    'molecular_crystals' : [1, 6, 7, 8],
    'silicon_bulk' : [14],
    'methane_liquid' : [1, 6],
    'methane_sulfonic' : [1,6,8,16],
}

gap_fit_params_default = {
    'cutoff': 5.0,
    'cutoff_transition_width': 1.0,
    #'n_max': 10,
    #'l_max': 12,
    'atom_width': 0.4,
    'soap_zeta': 2,
    'energy_scale': 1.0,
    'energy_reg': 0.001,
    'force_reg': 0.01
}

nl_sets = [
    {'n_max': 10, 'l_max': 12},
    {'n_max': 8, 'l_max': 6}
]

do_train_with_forces = [True, False]

gap_fit_params_fixed = {name: gap_fit_params_default.copy()
                        for name in system_names}

# Exceptions to the default params (why...?)
gap_fit_params_fixed['silicon_bulk']['atom_width'] = 0.5
gap_fit_params_fixed['silicon_bulk']['cutoff_transition_width'] = 1.0
gap_fit_params_fixed['silicon_bulk']['soap_zeta'] = 4
gap_fit_params_fixed['qm9']['atom_width'] = 0.3
gap_fit_params_fixed['qm9']['n_max'] = 12
gap_fit_params_fixed['qm9']['l_max'] = 9
gap_fit_params_fixed['molecular_crystals']['n_max'] = 9
gap_fit_params_fixed['molecular_crystals']['l_max'] = 9
gap_fit_params_fixed['methane_sulfonic']['n_max'] = 8
gap_fit_params_fixed['methane_sulfonic']['l_max'] = 6

n_sparse_all = [100, 200, 500, 1000, 2000, 5000, 9000]

for (system_name, param_set), n_sparse, nl_set, use_forces in itertools.product(
        gap_fit_params_fixed.items(), n_sparse_all, nl_sets, do_train_with_forces):
    param_set['system_name'] = system_name
    param_set['n_sparse'] = n_sparse
    # Shouldn't affect model evaluation speed, but check anyway
    # (for the two systems below)
    param_set['train_with_forces'] = use_forces
    # Only use the two sets of (n_max, l_max) parameters for these two systems
    if ((system_name == 'methane_liquid')
            or (system_name == 'silicon_bulk')):
        param_set['n_max'] = nl_set['n_max']
        param_set['l_max'] = nl_set['l_max']
    else:
        if param_set['train_with_forces']:
            continue
    job = project.open_job(param_set)
    job.init()
    # Metadata that is system-specific, but doesn't define a state point
    job.doc['system_sourcefile'] = system_filenames[system_name]
    job.doc['atoms_filename'] = os.path.basename(system_filenames[system_name])
    job.doc['energy_key'] = system_energy_keys[system_name]
    if use_forces:
        job.doc['force_key'] = system_force_keys[system_name]
    else:
        job.doc['force_key'] = 'none'
    job.doc['global_species'] = global_species[system_name]

