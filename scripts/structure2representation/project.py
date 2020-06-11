# project.py
from flow import FlowProject
from subprocess import Popen, PIPE
from os.path import join,dirname
from copy import deepcopy

import sys, os
import psutil


sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH
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

@FlowProject.operation
@FlowProject.post(ri_cpp_computed)
def compute_ri_cpp(job):
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
    # print([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])])
    p = Popen([group['executable'], job.fn(group['fn_in']), job.fn(group['fn_out'])], stderr=PIPE, stdout=PIPE)
    # print(p.stderr.read())
    data = fromjson(job.fn(group['fn_out']))
    data = data['results']
    tojson(job.fn(group['fn_res']), data)

@FlowProject.operation
@FlowProject.pre.after(compute_ri_cpp)
@FlowProject.post(lambda job: group['name'] in job.document)
def store_ri_cpp_in_document(job):
    data = fromjson(job.fn(group['fn_res']))
    job.document[group['name']] = data



if __name__ == '__main__':
    FlowProject().main()