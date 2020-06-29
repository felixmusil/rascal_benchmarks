# project.py
from flow import FlowProject
from subprocess import Popen, PIPE, run
from os.path import join,dirname
from copy import deepcopy
from ase.io import read
import sys, os
import numpy as np
from codetiming import Timer

sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
from utils.io import tojson, fromjson, fromfile, _decode
sys.path.insert(0, RASCAL_BUILD_PATH)
from rascal.representations import SphericalInvariants, SphericalExpansion
from rascal.representations.spherical_invariants import get_power_spectrum_index_mapping
from rascal.representations import SphericalInvariants
from rascal.models import Kernel, train_gap_model, SparsePoints, KRR
from rascal.neighbourlist import AtomsList
from rascal.utils import from_dict, CURFilter, fps, to_dict
from rascal.utils.random_filter import RandomFilter
from rascal.utils.io import dump_obj,load_obj

group = {
    'sparse_point_fn' : 'sparse_point.json',
    'knm_fn' : 'knm.npy',
    'kernel_fn' : 'kernel.json',
    'model_fn' : 'model.json',
    'benchmark_fn' : 'benchmark.json',
}

def get_feature_index_mapping(soap, species):
    u_species = np.unique(species)
    sp_pairs = soap.get_keys(u_species)
    n_max = soap.hypers['max_radial']
    l_max = soap.hypers['max_angular'] + 1
    if soap.hypers['soap_type'] == 'PowerSpectrum':
        feature_index_mapping = get_power_spectrum_index_mapping(
            sp_pairs, n_max, l_max)
    else:
        raise ValueError('Works only for PowerSpectrum')
    return feature_index_mapping

def get_randomly_sparsified_soap(data, rep):
    rep_hypers = deepcopy(data['representation'])
    feat_frac = rep_hypers['coefficient_subselection']
    global_species = rep_hypers['global_species']
    feature_mapping = get_feature_index_mapping(rep, global_species)
    ids = np.arange(len(feature_mapping))
    np.random.seed(data['seed'])
    np.random.shuffle(ids)
    sel_feat = {'a': [],'b': [], 'n1': [], 'n2': [], 'l': []}
    n_feat = int(feat_frac*len(feature_mapping))
    for idx in ids[:n_feat]:
        feat = feature_mapping[idx]
        for k,v in feat.items():
            sel_feat[k].append(int(v))
    rep_hypers['coefficient_subselection'] = sel_feat
    soap = SphericalInvariants(**rep_hypers)
    return soap,n_feat

@FlowProject.label
def sparse_point_computed(job):
    return job.isfile(group['sparse_point_fn'])
@FlowProject.label
def knm_computed(job):
    return job.isfile(group['knm_fn']) and job.isfile(group['kernel_fn'])
@FlowProject.label
def model_computed(job):
    return job.isfile(group['model_fn'])
@FlowProject.label
def benchmark_computed(job):
    return job.isfile(group['benchmark_fn'])


@FlowProject.operation
@FlowProject.post(sparse_point_computed)
def compute_sparse_point(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]
    soap = SphericalInvariants(**sp['representation'])
    managers = soap.transform(frames)
    compressor = RandomFilter(soap, **sp['RandomFilter'])
    X_pseudo = compressor.select_and_filter(managers)
    dump_obj(job.fn(group['sparse_point_fn']), X_pseudo)

def compute(i_frame,frame, soap, kernel, X_pseudo, grad=False):
    feat = soap.transform([frame])
    en_row = kernel(feat, X_pseudo)
    grad_rows = kernel(feat, X_pseudo, grad=(grad, False))
    return en_row, grad_rows

@FlowProject.operation
@FlowProject.pre.after(compute_sparse_point)
@FlowProject.post(knm_computed)
def compute_knm(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    X_pseudo = load_obj(job.fn(group['sparse_point_fn']))

    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = job.sp.train_with_grad
    soap = SphericalInvariants(**hypers)
    kernel = Kernel(soap, **sp['kernel'])

    Nstructures = len(frames)
    Ngrad_stride = [0]
    Ngrads = 0
    for frame in frames:
        n_at = len(frame)
        Ngrad_stride.append(n_at*3)
        Ngrads += n_at*3
    Ngrad_stride = np.cumsum(Ngrad_stride) + Nstructures
    dump_obj(job.fn(group['kernel_fn']), kernel)

    if job.sp.train_with_grad:
        KNM = np.zeros((Nstructures+Ngrads, X_pseudo.size()))
    else:
        KNM = np.zeros((Nstructures, X_pseudo.size()))


    for i_frame,frame in enumerate(frames):
        en_row, grad_rows = compute(i_frame,frame, soap, kernel, X_pseudo, grad=job.sp.train_with_grad)
        KNM[i_frame] = en_row
        if job.sp.train_with_grad:
            KNM[Ngrad_stride[i_frame]:Ngrad_stride[i_frame+1]] = grad_rows

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

def train_gap_model(kernel, frames, KNM_, X_pseudo, y_train, self_contributions, grad_train=None,
                    lambdas=None, jitter=1e-8):
    KMM = kernel(X_pseudo)
    Y = y_train.reshape((-1, 1)).copy()
    KNM = KNM_.copy()
    n_centers = Y.shape[0]
    Natoms = np.zeros(n_centers)
    Y0 = np.zeros((n_centers, 1))
    for iframe, frame in enumerate(frames):
        Natoms[iframe] = len(frame)
        numbers = frame.get_atomic_numbers()
        for sp in numbers:
            Y0[iframe] += self_contributions[sp]
    Y = Y - Y0
    delta = np.std(Y)
    print(delta)
    # lambdas[0] is provided per atom hence the '* np.sqrt(Natoms)'
    # the first n_centers rows of KNM are expected to refer to the
    # property
    KNM[:n_centers] /= lambdas[0] / delta * np.sqrt(Natoms)[:, None]
    Y /= lambdas[0] / delta * np.sqrt(Natoms)[:, None]

    if grad_train is not None:
        KNM[n_centers:] /= lambdas[1] / delta
        F = grad_train.reshape((-1, 1)).copy()
        F /= lambdas[1] / delta
        Y = np.vstack([Y, F])

    K = KMM + np.dot(KNM.T, KNM)
    eig,_ = np.linalg.eig(K)
    if eig.min() < 0:
        jitter = 1.05*abs(eig.min())
    else:
        jitter = 1e-9
    K[np.diag_indices_from(K)] += jitter

    Y = np.dot(KNM.T, Y)
    weights = np.linalg.solve(K, Y)
    model = KRR(weights, kernel, X_pseudo, self_contributions)

    # avoid memory clogging
    del K, KMM
    K, KMM = [], []

    return model

@FlowProject.operation
@FlowProject.pre.after(compute_knm)
@FlowProject.post(model_computed)
def train_model(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    X_pseudo = load_obj(job.fn(group['sparse_point_fn']))
    kernel = load_obj(job.fn(group['kernel_fn']))
    KNM = np.load(job.fn(group['knm_fn']), mmap_mode='r')

    y,f = extract_ref(frames, info_key='dft_energy',array_key='zeros')
    if not job.sp.train_with_grad:
        f = None
    else:
        f = -f

    model = train_gap_model(kernel, frames, KNM, X_pseudo, y, sp['self_contributions'],
                            grad_train=f, lambdas=[1e-3], jitter=1e-8)

    dump_obj(job.fn(group['model_fn']), model)


@FlowProject.operation
@FlowProject.pre.after(train_model)
@FlowProject.post(benchmark_computed)
def compute_benchmark(job):
    sp = _decode(job.statepoint())
    st,lg = job.sp.start_structure, job.sp.n_structures
    frames = fromfile(job.sp.filename)[st:st+lg]

    model = load_obj(job.fn(group['model_fn']))
    soap = model.get_representation_calculator()

    hypers = soap._get_init_params()
    hypers['compute_gradients'] = True
    soap = SphericalInvariants(**hypers)

    rc = sp['representation']['interaction_cutoff']
    nl_options = [
            dict(name='centers', args=[]),
            dict(name='neighbourlist', args=dict(cutoff=rc)),
            # dict(name='halflist', args=dict()),
            dict(name="centercontribution", args=dict()),
            dict(name='strict', args=dict(cutoff=rc))
        ]

    kernel = Kernel(soap, **sp['kernel'])

    N_ITERATIONS = sp['N_ITERATIONS']
    N_ITERATIONS = 2

    t1 = Timer(name='NL',text='')
    t2 = Timer(name='rep with grad',text='')
    t2b = Timer(name='rep no grad',text='')
    t3 = Timer(name='kernel',text='')
    t3b = Timer(name='kernel grad',text='')
    t4 = Timer(name='pred energy',text='')
    t5 = Timer(name='pred forces',text='')

    n_atoms = 0
    for frame in frames:
        n_atoms += len(frame)

    for _ in range(N_ITERATIONS):
        with t1:
            managers = AtomsList(frames, nl_options)
        with t2:
            managers = soap.transform(managers)

        Y0 = model._get_property_baseline(managers)

        with t3:
            KNM = kernel(managers, model.X_train, (False, False))
        with t4:
            Y0 + np.dot(KNM, model.weights).reshape((-1))

        with t3b:
            KNM = kernel(managers, model.X_train, (True, False))
        with t5:
            np.dot(KNM, model.weights).reshape((-1, 3))

    print(Timer.timers)
    print(Timer.timers.mean('rep with grad'),Timer.timers.stdev('rep with grad'))

# @FlowProject.operation
# @FlowProject.pre.after(compute_si_cpp)
# @FlowProject.post(lambda job: group['name'] in job.document)
# def store_si_cpp_in_document(job):
#     data = fromjson(job.fn(group['fn_res']))
#     job.document = data


if __name__ == '__main__':
    FlowProject().main()