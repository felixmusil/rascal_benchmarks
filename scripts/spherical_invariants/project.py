# project.py
from flow import FlowProject
from subprocess import Popen, PIPE
from os.path import join,dirname
from copy import deepcopy
from memory_profiler import memory_usage
from ase.io import read
import sys, os
import numpy as np
from time import ctime

sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
from utils.io import tojson, fromjson
sys.path.insert(0, RASCAL_BUILD_PATH)
from rascal.representations import SphericalInvariants, SphericalExpansion
from rascal.representations.spherical_invariants import get_power_spectrum_index_mapping

groups = {
    'spherical_invariants_cpp' : {
        'name' : 'spherical_invariants_cpp',
        'fn_out' : 'out_si_cpp.json',
        'fn_res' : 'res_si_cpp.json',
        'fn_in' : 'in_si_cpp.json',
        'executable' : {
            'Full':join(BUILD_PATH,'src/spherical_invariants'),
            'Half':join(BUILD_PATH,'src/spherical_invariants_half'),
            },
    },
}

group = groups['spherical_invariants_cpp']

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
def si_cpp_computed(job):
    return job.isfile(group['fn_out']) and job.isfile(group['fn_res'])

@FlowProject.operation
@FlowProject.post(si_cpp_computed)
def compute_si_cpp(job):
    # setup input for the script
    data = job.statepoint()
    np.random.seed(data['seed'])

    if data['representation']['coefficient_subselection'] is not None:
        rpr = deepcopy(data['representation'])
        rpr.pop('coefficient_subselection')
        rep = SphericalInvariants(**rpr)
        rep,n_feat = get_randomly_sparsified_soap(data, rep)
    else:
        n_feat = None
        rep = SphericalInvariants(**data['representation'])

    data['calculator'] = rep.hypers
    tojson(job.fn(group['fn_in']), data)
    # look at memory footprint
    p = Popen([group['executable'][data['nl_type']], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    max_mem = memory_usage(p, interval=0.1, max_usage=True)
    # look at timings
    p = Popen([group['executable'][data['nl_type']], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    if p.stderr.read(): print(p.stderr.read())
    data = fromjson(job.fn(group['fn_out']))
    data = data['results']
    data['n_features'] = n_feat
    data['mem_max'] = max_mem
    data['mem_unit'] = 'MiB'
    tojson(job.fn(group['fn_res']), data)

@FlowProject.operation
@FlowProject.pre.after(compute_si_cpp)
@FlowProject.post(lambda job: group['name'] in job.document)
def store_si_cpp_in_document(job):
    data = fromjson(job.fn(group['fn_res']))
    job.document = data


if __name__ == '__main__':
    print(ctime())
    FlowProject().main()
    print(ctime())