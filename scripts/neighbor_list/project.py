# project.py
from flow import FlowProject
from subprocess import Popen, PIPE
from os.path import join,dirname
from copy import deepcopy
from memory_profiler import memory_usage


import sys, os


sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
from utils.io import tojson, fromjson
# sys.path.insert(0, RASCAL_BUILD_PATH)
# from rascal.representations import SphericalInvariants, SphericalExpansion


groups = {
    'neighbor_list_cpp' : {
        'name' : 'neighbor_list_cpp',
        'fn_out' : 'out_nl_cpp.json',
        'fn_res' : 'res_nl_cpp.json',
        'fn_in' : 'in_nl_cpp.json',
        'executable' : join(BUILD_PATH,'src/neighbor_list'),
    },
}

group = groups['neighbor_list_cpp']

@FlowProject.label
def nl_cpp_computed(job):
    return job.isfile(group['fn_out']) and job.isfile(group['fn_res'])

@FlowProject.operation
@FlowProject.post(nl_cpp_computed)
def compute_nl_cpp(job):
    # setup input for the script
    data = job.statepoint()

    tojson(job.fn(group['fn_in']), data)
    # look at timings
    p = Popen([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])], stdout=PIPE, stderr=PIPE)
    if p.stderr.read(): print(p.stderr.read())
    data = fromjson(job.fn(group['fn_out']))
    data = data['results']
    tojson(job.fn(group['fn_res']), data)

@FlowProject.operation
@FlowProject.pre.after(compute_nl_cpp)
@FlowProject.post(lambda job: group['name'] in job.document)
def store_nl_cpp_in_document(job):
    data = fromjson(job.fn(group['fn_res']))
    job.document = data

if __name__ == '__main__':
    FlowProject().main()