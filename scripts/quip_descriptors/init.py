"""Initialize the data space for QUIP overall model benchmarks"""

import itertools
import os

import signac

project = signac.init_project('quip-descriptors')

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

global_species = {
    'qm9' : [1,6,7,8,9],
    'molecular_crystals' : [1, 6, 7, 8],
    'silicon_bulk' : [14],
    'methane_liquid' : [1, 6],
    'methane_sulfonic' : [1,6,8,16],
}

test_subsets = {
    'qm9': (0, 800),
    'molecular_crystals': (0, 90),
    'silicon_bulk': (600, 120),
    'methane_liquid': (100, 100),
    'methane_sulfonic': (0, 50)
}

soap_params_default = {
    'cutoff': 5.0,
    'cutoff_transition_width': 1.0,
    'n_max': 10,
    'l_max': 12,
    'atom_width': 0.4,
    'soap_zeta': 2,
}

# Exceptions to the default params (why...?)
soap_params_fixed = {name: soap_params_default.copy()
                        for name in system_names}
soap_params_fixed['silicon_bulk']['atom_width'] = 0.5
soap_params_fixed['silicon_bulk']['cutoff_transition_width'] = 1.0
soap_params_fixed['silicon_bulk']['soap_zeta'] = 4
soap_params_fixed['qm9']['atom_width'] = 0.3
soap_params_fixed['qm9']['n_max'] = 12
soap_params_fixed['qm9']['l_max'] = 9
soap_params_fixed['molecular_crystals']['n_max'] = 9
soap_params_fixed['molecular_crystals']['l_max'] = 9
soap_params_fixed['methane_sulfonic']['n_max'] = 8
soap_params_fixed['methane_sulfonic']['l_max'] = 6

n_max_sweep = range(2,16,2)
l_max_sweep = range(1,16,2)

def populate_job(job, system_name):
    """Fill out metadata that is system-specific, but doesn't define a state point"""
    job.doc['system_sourcefile'] = system_filenames[system_name]
    job.doc['atoms_filename'] = os.path.basename(system_filenames[system_name])
    job.doc['global_species'] = global_species[system_name]
    job.doc['test_subset'] = test_subsets[system_name]

# First with the default (n, l) params
for system_name, param_set in soap_params_fixed.items():
    param_set['system_name'] = system_name
    job = project.open_job(param_set)
    job.init()
    populate_job(job, system_name)

# Then with a whole Cartesian (n, l) parameter space sweep
for (system_name, param_set), n_max, l_max in itertools.product(
        soap_params_fixed.items(), n_max_sweep, l_max_sweep):
    param_set['system_name'] = system_name
    param_set['n_max'] = n_max
    param_set['l_max'] = l_max
    job = project.open_job(param_set)
    job.init()
    populate_job(job, system_name)
