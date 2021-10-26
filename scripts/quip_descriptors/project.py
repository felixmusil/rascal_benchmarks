"""Run SOAP benchmarks using QUIP"""

import errno
import logging
import os
import re
import subprocess

import numpy as np
import signac
import flow
from flow import FlowProject, directives

# Sigh... the relative-import machinery just doesn't work with signac's
# project layout, does it.
#from ..utils.quip import store_timing_data
import sys
sys.path.append(os.path.join(os.getcwd(), '..'))
from utils.quip import store_timing_data


LOGGER = logging.getLogger(__name__)


class HelvetiosEnvironment(flow.environment.DefaultSlurmEnvironment):

    hostname_pattern = r'^helvetios$'
    template = 'helvetios.sh'

    @classmethod
    def add_args(cls, parser):
        super(HelvetiosEnvironment, cls).add_args(parser)
        parser.add_argument('--job-output',
                            help="Write output to the specified file. "
                            "Slurm substitutions (e.g. %%j) are accepted.")
        parser.add_argument('--partition',
                            help="Name of the slurm partition to submit to")


@FlowProject.label
def atoms_file_available(job):
    return job.isfile(os.path.join(signac.get_project().root_directory(),
                                   'xyz_files',
                                   job.doc.atoms_filename))


@FlowProject.label
def has_soap_timings(job):
    return ('energy_ip_total_time_peratom_mean' in job.doc)


@FlowProject.label
def has_soap_grad_timings(job):
    return ('force_ip_total_time_peratom_mean' in job.doc)


def has_libatoms_exit_message(filename):
    """Check if the (potentially very large) file ends in the expected string"""
    fileend = "libAtoms::Finalise: Bye-Bye!\n"
    # Warning: Potential race condition(s)
    # (also implicitly assuming UTF-8 encoding)
    fsize = os.stat(filename).st_size
    if fsize < len(fileend):
        return False
    with open(filename, 'r') as outfile:
        outfile.seek(fsize - len(fileend))
        return (outfile.read(len(fileend)) == fileend)


def has_soap_output(job):
    outfile_name = 'soap_output.out'
    if not job.isfile(outfile_name):
        return False
    else:
        return has_libatoms_exit_message(job.fn(outfile_name))


def has_soap_grad_output(job):
    outfile_name = 'soap_grad_output.out'
    if not job.isfile(outfile_name):
        return False
    else:
        return has_libatoms_exit_message(job.fn(outfile_name))


def build_quip_command_line(job, do_gradients=False):
    cmd = 'quip'
    args = []
    # May use a separate test atoms file later
    args.append('atoms_filename={:s}'.format(
        os.path.join(signac.get_project().root_directory(),
                     'xyz_files',
                     job.doc.atoms_filename)))
    args.append(
        'descriptor_str={{soap atom_sigma={sp.atom_width:f} l_max={sp.l_max:d} '
        'n_max={sp.n_max:d} cutoff={sp.cutoff:f} '
        'cutoff_transition_width={sp.cutoff_transition_width:f} '
        'energy_scale={sp.energy_scale:f} add_species n_species={n_species:d} '
        'species_Z={{{{ {species_str:s} }}}} n_sparse={sp.n_sparse:d} '
        'sparse_method=cur_points covariance_type=dot_product '
        'soap_exponent={sp.soap_zeta:d} }}'.format(
            sp=job.sp, n_species=len(job.doc.global_species),
            species_str=(' '.join(str(num) for num in job.doc.global_species))))
    if do_gradients:
        args.append('do_grad_descriptor')
    args.append('timing')
    return [cmd, ] + args


timings = FlowProject.make_group(name='timings')


@timings
@FlowProject.operation
@FlowProject.pre(atoms_file_available)
@FlowProject.post(has_soap_output)
@directives(omp_num_threads=1)
def time_quip_soap(job):
    """Run QUIP descriptor timings"""
    old_dir = os.getcwd()
    quip_cmdline = build_quip_command_line(job)
    print(quip_cmdline)
    try:
        os.chdir(job.workspace())
        with open('soap_output.out', 'w') as outfile:
            subprocess.run(quip_cmdline, stdout=outfile)
    finally:
        os.chdir(old_dir)


@timings
@FlowProject.operation
@FlowProject.pre(atoms_file_available)
@FlowProject.post(has_soap_grad_output)
@directives(omp_num_threads=1)
def time_quip_soap_gradients(job):
    """Run QUIP descriptor timings with gradients"""
    old_dir = os.getcwd()
    quip_cmdline = build_quip_command_line(job, do_forces=True)
    print(quip_cmdline)
    try:
        os.chdir(job.workspace())
        with open('soap_grad_output.out', 'w') as outfile:
            subprocess.run(quip_cmdline, stdout=outfile)
    finally:
        os.chdir(old_dir)


#TODO merge with quip-potential version from quip-benchmarks project
def parse_quip_timings(job, filename):
    """Parse timing data from 'quip' output and check for consistency"""
    with open(job.fn(filename), 'r') as outfile:
        connect_timer = []
        soap_timer = []
        gap_timer = []
        ip_timer = []
        natoms_list = []
        for line in outfile:
            atmatch = re.match(r'AT (\d+)', line)
            if atmatch:
                natoms_list.append(int(atmatch.group(1)))
            if re.match('TIMER: soap_calc', line):
                soap_timer.append(float(line.split()[7]))
    # Check for consistency
    n_configs = len(natoms_list)
    # Check index file: Not fatal if wrong or missing, but should be noticed
    xyz_index_file = os.path.join(signac.get_project().root_directory(),
                                  'xyz_files',
                                  job.doc.atoms_filename + '.idx')
    if os.path.isfile(xyz_index_file):
        with open(xyz_index_file, 'r') as idx_f:
            n_configs_index = int(idx_f.readline())
        if n_configs_index != n_configs:
            LOGGER.warning("Mismatch between {:d} configurations in energy "
                           "output and {:d} in index file {:s}".format(
                               n_configs, n_configs_index, xyz_index_file))
    else:
        LOGGER.warning("Could not find index file {:s}".format(xyz_index_file))
    error_message = ("Wrong number {number:d} of {time_type:s} timings in file"
                     " {filename:s}")
    timers = (connect_timer, soap_timer, gap_timer, ip_timer)
    labels = ('calc_connect', 'soap', 'gap', 'ip_total')
    n_soaps = len(job.doc['global_species'])
    multiplicities = (2, n_soaps, n_soaps, 1)
    timings_out = {'natoms': np.array(natoms_list)}
    if len(soap_timer) != n_soaps*n_configs:
        raise ValueError(error_message.format(
            number=len(soap_timer),
            time_type='soap',
            filename=job.fn(filename)
        ))
    else:
        timings_out['soap'] = np.array(soap_timer).reshape((-1, n_soaps))
    return timings_out


@FlowProject.operation
@FlowProject.pre(has_soap_output)
@FlowProject.post(has_soap_timings)
@FlowProject.post.never
def process_energy_timings(job):
    timing_data = parse_quip_timings(job, 'soap_output.out')
    store_timing_data(job, timing_data, 'soap')


@FlowProject.operation
@FlowProject.pre(has_soap_grad_output)
@FlowProject.post(has_soap_grad_timings)
@FlowProject.post.never
def process_force_timings(job):
    timing_data = parse_quip_timings(job, 'soap_grad_output.out')
    store_timing_data(job, timing_data, 'soap_grad')


if __name__ == "__main__":
    logging.basicConfig()
    FlowProject().main()
