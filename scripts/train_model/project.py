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
from concurrent.futures import as_completed, ProcessPoolExecutor, ThreadPoolExecutor
from sklearn.model_selection._split import KFold

group = {
    'sparse_point_fn' : 'sparse_point.json',
    'feature_fn' : 'feature.json',
    'knm_fn' : 'knm.pck',
    'kernel_fn' : 'kernel.json',
    'model_fn' : 'model.json',
    'cv_score': 'cv_score.json',
    'fps_features' : '../../feature_selection/{}/fps_features.pck',
    'fps_samples' : '../../feature_selection/{}/fps_samples.pck',
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
    'methane_liquid': 'force',
}




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
@FlowProject.post(sparse_point_computed)
def compute_sparse_point(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]
    soap = SphericalInvariants(**sp['representation'])
    managers = soap.transform(frames)

    fcompressor = frompickle(group['fps_features'].format(sp['name']))
    feature_subselection = fcompressor.filter(managers, sp['feature_subselection'])

    sp['representation']['coefficient_subselection'] = feature_subselection['coefficient_subselection']
    soap = SphericalInvariants(**sp['representation'])
    managers = soap.transform(frames)

    compressor = frompickle(group['fps_samples'].format(sp['name']))
    compressor._representation = soap
    X_pseudo = compressor.filter(managers, sp['sparse_point_subselection'])
    dump_obj(job.fn(group['sparse_point_fn']), X_pseudo)

def init_proc(sp_fn):
    global X_pseudos
    X_pseudos = load_obj(sp_fn)


def compute(i_frame,frame,zeta,train_with_grad):
    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = train_with_grad
    soap = SphericalInvariants(**hypers)
    kernel = Kernel(soap, name='GAP', zeta=zeta, target_type='Structure', kernel_type='Sparse')
    feat = soap.transform([frame])
    en_row = kernel(feat, X_pseudo)
    grad_rows = None
    if train_with_grad:
        grad_rows = kernel(feat, X_pseudo, grad=(True, False))
    return en_row, grad_rows

@FlowProject.operation
@FlowProject.pre.after(compute_sparse_point)
@FlowProject.post(knm_computed)
def compute_knm(job):

    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    Nstructures = len(frames)
    Ngrad_stride = [0]
    Ngrads = 0
    for frame in frames:
        n_at = len(frame)
        Ngrad_stride.append(n_at*3)
        Ngrads += n_at*3
    Ngrad_stride = np.cumsum(Ngrad_stride) + Nstructures

    X_pseudo = load_obj(job.fn(group['sparse_point_fn']))
    sp_fn = job.fn(group['sparse_point_fn'])
    n_sparse = X_pseudo.size()

    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = job.sp.train_with_grad
    zeta = job.sp['kernel']['zeta']

    with ProcessPoolExecutor(max_workers=OMP_NUM_THREADS, initializer=init_proc, initargs=(sp_fn)) as executor:
        if job.sp.train_with_grad:
            KNM = np.zeros((Nstructures+Ngrads, n_sparse))
        else:
            KNM = np.zeros((Nstructures, n_sparse))

        future_to_compute = {executor.submit(compute, i_frame, frame, zeta, job.sp.train_with_grad):i_frame
                                                     for i_frame,frame in enumerate(frames)}
        pbar = tqdm(total=len(future_to_compute),desc='KNM')
        for future in as_completed(future_to_compute):
            i_frame = future_to_compute[future]
            en_row, grad_rows = future.result()
            KNM[i_frame] = en_row
            if job.sp.train_with_grad:
                KNM[Ngrad_stride[i_frame]:Ngrad_stride[i_frame+1]] = grad_rows
            pbar.update()
        np.save(job.fn(group['knm_fn']), KNM)


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
def cv_score_model(job):
    cv = KFold(n_splits=4, shuffle=True, random_state=10)
    lamda_es = [1e-1,5e-2,1e-2,5e-3,1e-3]
    lamda_fs = [0.5,1e-1,5e-2,1e-2,5e-3,1e-3]

    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    X_pseudo = load_obj(job.fn(group['sparse_point_fn']))

    n_sparse = X_pseudo.size()
    n_feature = sp['sparse_point_subselection']

    KNM = np.load(job.fn(group['knm_fn']), mmap_mode=None)
    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = job.sp.train_with_grad
    soap = SphericalInvariants(**hypers)
    kernel = Kernel(soap, **job.sp['kernel'])

    scores_ = get_cv_scores(cv, KNM, frames, kernel, X_pseudo,
                           sp['self_contributions'], lamda_es, lamda_fs, **{'n_sparse':n_sparse, 'n_feature':n_feature})
    tojson(group['cv_score'], scores_)


@FlowProject.operation
@FlowProject.pre.after(train_model_)
@FlowProject.post(lambda job: 'done' in job.document)
def store_benchmark_in_document(job):
    job.document = {'done':True}


if __name__ == '__main__':
    print(ctime())
    FlowProject().main()
    print(ctime())