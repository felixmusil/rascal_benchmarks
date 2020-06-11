# init.py
import signac
import sys
from os.path import join,dirname
print(dirname(__file__))
sys.path.insert(0, join(dirname(__file__), '../'))
from representation_params import REPRESENTATION_PARAMS

project = signac.init_project('structure2representation')


for REPRESENTATION_PARAM in REPRESENTATION_PARAMS:
    job = project.open_job(REPRESENTATION_PARAM)
    job.init()