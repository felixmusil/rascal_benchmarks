# project.py
from flow import FlowProject
from subprocess import Popen, PIPE
from os.path import join,dirname
from copy import deepcopy
from memory_profiler import memory_usage
from time import ctime

import sys, os

sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
from utils.io import tojson, fromjson
sys.path.insert(0, RASCAL_BUILD_PATH)
from rascal.representations import SphericalInvariants, SphericalExpansion


groups = {
    'spherical_expansion_cpp' : {
        'name' : 'spherical_expansion_cpp',
        'fn_out' : 'out_se_cpp.json',
        'fn_res' : 'res_se_cpp.json',
        'fn_in' : 'in_se_cpp.json',
        'executable' : {
            'Full':join(BUILD_PATH,'src/spherical_expansion'),
            'Half':join(BUILD_PATH,'src/spherical_expansion_half'),
            },
    },
}

group = groups['spherical_expansion_cpp']

@FlowProject.label
def se_cpp_computed(job):
    return job.isfile(group['fn_out']) and job.isfile(group['fn_res'])

@FlowProject.operation
@FlowProject.post(se_cpp_computed)
def compute_se_cpp(job):
    # setup input for the script
    data = job.statepoint()
    rep = SphericalExpansion(**data['representation'])
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
    data['mem_max'] = max_mem
    data['mem_unit'] = 'MiB'
    tojson(job.fn(group['fn_res']), data)

@FlowProject.operation
@FlowProject.pre.after(compute_se_cpp)
@FlowProject.post(lambda job: group['name'] in job.document)
def store_se_cpp_in_document(job):
    data = fromjson(job.fn(group['fn_res']))
    job.document = data


if __name__ == '__main__':
    print(ctime())
    FlowProject().main()
    print(ctime())