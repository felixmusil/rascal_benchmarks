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
    'radial_integral_cpp' : {
        'name' : 'radial_integral_cpp',
        'fn_out' : 'out_ri_cpp.json',
        'fn_res' : 'res_ri_cpp.json',
        'fn_in' : 'in_ri_cpp.json',
        'executable' : join(BUILD_PATH,'src/radial_integral'),
    },
}

group = groups['radial_integral_cpp']

@FlowProject.label
def ri_cpp_computed(job):
    return job.isfile(group['fn_out']) and job.isfile(group['fn_res'])

def run(job, N_ITERATION, max_mem):
    data = job.statepoint()
    data['N_ITERATIONS'] = N_ITERATION
    tojson(job.fn(group['fn_in']), data)
    p = Popen([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    # if p.stderr.read(): print(p.stderr.read())
    data = fromjson(job.fn(group['fn_out']))
    data = data['results']
    data['mem_max'] = max_mem
    data['mem_unit'] = 'MiB'
    tojson(job.fn(group['fn_res']), data)
    # job.document.elapsed_mean < 3*job.document.elapsed_std
    return data['elapsed_mean'] > 4*data['elapsed_std']


@FlowProject.operation
@FlowProject.post(ri_cpp_computed)
def compute_ri_cpp(job):
    # setup input for the script
    data = job.statepoint()
    rep = SphericalExpansion(**data['representation'])
    data['calculator'] = rep.hypers

    tojson(job.fn(group['fn_in']), data)
    # look at memory footprint
    p = Popen([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    max_mem = memory_usage(p, interval=0.001, max_usage=True)
    # look at timings
    carry_on = True
    N_ITERATIONS = [int(data['N_ITERATIONS']*(1+f)) for f in [0., 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100]]
    for N_ITERATION in N_ITERATIONS:
        if run(job, N_ITERATION, max_mem):
            break


@FlowProject.operation
@FlowProject.pre.after(compute_ri_cpp)
@FlowProject.post(lambda job: group['name'] in job.document)
def store_ri_cpp_in_document(job):
    data = fromjson(job.fn(group['fn_res']))
    job.document = data


if __name__ == '__main__':
    print(ctime())
    FlowProject().main()
    print(ctime())