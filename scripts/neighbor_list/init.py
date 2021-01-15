# init.py
import signac
import sys
from copy import deepcopy
import numpy as np
from os.path import join,dirname
from itertools import product

sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH
# from representation_params import REPRESENTATION_PARAMS

project = signac.init_project('neighbor_list')

rcs = [3,4,5,6,7,8]
names = ['qm9', 'molecular_crystals', 'silicon_bulk']

fns = {
    'qm9': 'qm9.json',
    'molecular_crystals': 'molecular_crystals_50.json',
    'silicon_bulk': 'silicon_bulk.json',
    'methane_liquid': 'methane_liquid.json',
    'methane_sulfonic': 'methane_sulfonic_150.json'
}

misc_entries = {
    'qm9' : {
        "N_ITERATIONS" : 500,
        "n_structures" : 500,
        'start_structure': 0,
    },
    'molecular_crystals' : {
        "N_ITERATIONS" : 500,
        'start_structure': 0,
        "n_structures" : 50,
    },
    'silicon_bulk' : {
        "N_ITERATIONS" : 500,
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


for name, rc in product(names, rcs):
    rep_args = {
        "filename" : join(STRUCTURE_PATH,fns[name]),
        "adaptors" : [
            {"initialization_arguments": {"cutoff": rc}, "name": "neighbourlist"},
            {"initialization_arguments": {}, "name": "centercontribution"},
            {"initialization_arguments": {"cutoff": rc}, "name": "strict"}
        ],
    }
    rep_args.update(**misc_entries[name])

    job = project.open_job(rep_args)
    job.init()