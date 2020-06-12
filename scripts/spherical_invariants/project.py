# project.py
from flow import FlowProject
from subprocess import Popen, PIPE
from os.path import join,dirname
from copy import deepcopy
from memory_profiler import memory_usage

import sys, os

sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
from utils.io import tojson, fromjson
sys.path.insert(0, RASCAL_BUILD_PATH)
from rascal.representations import SphericalInvariants, SphericalExpansion


groups = {
    'spherical_invariants_cpp' : {
        'name' : 'spherical_invariants_cpp',
        'fn_out' : 'out_si_cpp.json',
        'fn_res' : 'res_si_cpp.json',
        'fn_in' : 'in_si_cpp.json',
        'executable' : join(BUILD_PATH,'src/spherical_invariants'),
    },
}

group = groups['spherical_invariants_cpp']

@FlowProject.label
def si_cpp_computed(job):
    return job.isfile(group['fn_out']) and job.isfile(group['fn_res'])

@FlowProject.operation
@FlowProject.post(si_cpp_computed)
def compute_si_cpp(job):
    # setup input for the script
    data = job.statepoint()
    rep = SphericalInvariants(**data['representation'])
    data['calculator'] = rep.hypers
    cutoff = rep.hypers['cutoff_function']['cutoff']['value']
    data['adaptors'] = [
        {"initialization_arguments": {"cutoff": cutoff}, "name":   "neighbourlist"},
        {"initialization_arguments": {}, "name": "centercontribution"},
        {"initialization_arguments": {"cutoff": cutoff}, "name": "strict"}
    ]
    tojson(job.fn(group['fn_in']), data)
    # look at memory footprint
    p = Popen([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    max_mem = memory_usage(p, interval=0.001, max_usage=True)
    # look at timings
    p = Popen([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    data = fromjson(job.fn(group['fn_out']))
    data = data['results']
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
    FlowProject().main()