import os
OMP_NUM_THREADS = 12
os.environ["OMP_NUM_THREADS"] = "{}".format(OMP_NUM_THREADS)

# project.py
from flow import FlowProject
from subprocess import Popen, PIPE, run
from os.path import join,dirname
from copy import deepcopy
from ase.io import read
import sys, os
import numpy as np
from time import sleep, ctime

sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
from utils.io import tojson, fromjson, fromfile, _decode, topickle, frompickle
from utils.model import KnmPool, train_gap_model
from utils import Timer
sys.path.insert(0, RASCAL_BUILD_PATH)
from rascal.representations import SphericalInvariants, SphericalExpansion
from rascal.representations.spherical_invariants import get_power_spectrum_index_mapping
from rascal.representations import SphericalInvariants
from rascal.models import Kernel, SparsePoints
from rascal.neighbourlist import AtomsList
from rascal.utils import from_dict, CURFilter, FPSFilter, fps, to_dict
from rascal.utils.random_filter import RandomFilter
from rascal.utils.io import dump_obj,load_obj
from ase.build import make_supercell
from tqdm import tqdm



group = {
    'sparse_point_fn' : 'sparse_point.json',
    'feature_fn' : 'feature.json',
    'knm_fn' : 'knm.pck',
    'kernel_fn' : 'kernel.json',
    'model_fn' : 'model.json',
}

en_field = {
    'qm9': 'dHf_peratom',
    'molecular_crystals': 'ENERGY',
    'silicon_bulk': 'dft_energy',
    'methane_liquid': 'energy',
    'methane_sulfonic': 'energy',
}

force_field = {
    'silicon_bulk': 'dft_force',
}



@FlowProject.label
def feature_sel_computed(job):
    return job.isfile(group['feature_fn'])
@FlowProject.label
def sparse_point_computed(job):
    return job.isfile(group['sparse_point_fn'])
@FlowProject.label
def knm_computed(job):
    return job.isfile(group['knm_fn'])

@FlowProject.label
def model_computed(job):
    return job.isfile(group['model_fn'])


@FlowProject.operation
@FlowProject.post(feature_sel_computed)
def compute_feature_selection(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]
    soap = SphericalInvariants(**sp['representation'])
    managers = soap.transform(frames)
    compressor = FPSFilter(soap, **sp['feature_subselection'])
    compressor.select(managers)
    feature_subselection = compressor.filter(managers)

    tojson(job.fn(group['feature_fn']), feature_subselection)

@FlowProject.operation
@FlowProject.pre.after(compute_feature_selection)
@FlowProject.post(sparse_point_computed)
def compute_sparse_point(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]
    soap = SphericalInvariants(**sp['representation'])
    managers = soap.transform(frames)
    compressor = FPSFilter(soap, **sp['sparse_point_subselection'])
    compressor.select(managers)

    feature_subselection = fromjson(job.fn(group['feature_fn']))
    sp['representation']['coefficient_subselection'] = feature_subselection['coefficient_subselection']
    soap = SphericalInvariants(**sp['representation'])
    managers = soap.transform(frames)
    compressor._representation = soap
    X_pseudo = compressor.filter(managers)
    dump_obj(job.fn(group['sparse_point_fn']), X_pseudo)


@FlowProject.operation
@FlowProject.pre.after(compute_sparse_point)
@FlowProject.post(knm_computed)
def compute_knm(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    X_pseudo = load_obj(job.fn(group['sparse_point_fn']))
    sp_fn = job.fn(group['sparse_point_fn'])

    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = job.sp.train_with_grad
    if job.sp.train_with_grad:
        pool = KnmPool(ncpu=OMP_NUM_THREADS, energy_tag=en_field[job.sp.name], forces_tag=force_field[job.sp.name])
    else:
        pool = KnmPool(ncpu=OMP_NUM_THREADS, energy_tag=en_field[job.sp.name])

    zeta = job.sp['kernel']['zeta']
    res = pool.run(frames, zeta=zeta, sparse_points_fn=sp_fn, self_contributions=sp['self_contributions'])

    topickle(job.fn(group['knm_fn']), res)

def extract_ref(frames,info_key='energy',array_key='forces'):
    y,f = [], []
    for frame in frames:
        y.append(frame.info[info_key])
        if array_key is None:
            pass
        elif array_key == 'zeros':
            f.append(np.zeros(frame.get_positions().shape))
        else:
            f.append(frame.get_array(array_key))
    y= np.array(y)
    try:
        f = np.concatenate(f)
    except:
        pass
    return y,f


@FlowProject.operation
@FlowProject.pre.after(compute_knm)
@FlowProject.post(model_computed)
def train_model_(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    X_pseudo = load_obj(job.fn(group['sparse_point_fn']))
    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = True
    soap = SphericalInvariants(**hypers)
    kernel = Kernel(soap, **job.sp['kernel'])

    res = frompickle(job.fn(group['knm_fn']))

    if job.sp.train_with_grad:
        model = train_gap_model(kernel, X_pseudo, res['energy'],
                                sp['self_contributions'], grads=res['grads'],
                                lambdas=[1e-3, 1e-2], jitter=1e-6)
    else:
        model = train_gap_model(kernel, X_pseudo, res['energy'],
                                sp['self_contributions'],
                                lambdas=[1e-3], jitter=1e-6)

    dump_obj(job.fn(group['model_fn']), model)


@FlowProject.operation
@FlowProject.pre.after(train_model_)
@FlowProject.post(lambda job: 'done' in job.document)
def store_benchmark_in_document(job):
    job.document = {'done':True}


if __name__ == '__main__':
    print(ctime())
    FlowProject().main()
    print(ctime())