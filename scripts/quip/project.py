"""Run QUIP benchmarks"""

import os
import errno
import logging
import subprocess

import numpy as np
import signac
import flow
from flow import FlowProject, directives


LOGGER = logging.getLogger(__name__)
LOGGER.basicConfig()


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


#TODO is this a workflow operation?
def build_gap_fit_command_line(job):
    cmd = 'gap_fit'
    args = []
    args.append('at_file={:s}'.format(
        os.path.join(signac.get_project().root_directory(),
                     'xyz_files',
                     job.doc.atoms_filename)))
    args.append(
        'gap={{soap atom_sigma={sp.atom_width:f} l_max={sp.l_max:d} '
        'n_max={sp.n_max:d} cutoff={sp.cutoff:f} '
        'cutoff_transition_width={sp.cutoff_transition_width:f} '
        'energy_scale={sp.energy_scale:f} add_species n_species={n_species:d} '
        'species_Z={{{{ {species_str:s} }}}} n_sparse={sp.n_sparse:d} '
        'sparse_method=cur_points covariance_type=dot_product '
        'soap_exponent={sp.soap_zeta:d} }}'.format(
            sp=job.sp, n_species=len(job.doc.global_species),
            species_str=(' '.join(str(num) for num in job.doc.global_species))))
    args.append(
        'default_kernel_regularisation={{ {0.energy_reg:f} {0.force_reg:f}'
        ' 1.0 1.0 }}'.format(job.sp))
    args.append('energy_parameter_name={0.energy_key:s}'.format(job.doc))
    args.append('force_parameter_name={0.force_key:s}'.format(job.doc))
    # We're not fitting with virials, no matter what
    args.append('virial_parameter_name=none')
    args.append('e0_method=average')
    args.append('gp_file={:s}'.format(job.fn('potential.xml')))
    return [cmd, ] + args


def build_quip_command_line(job, do_forces=False):
    cmd = 'quip'
    args = []
    # May use a separate test atoms file later
    args.append('atoms_filename={:s}'.format(
        os.path.join(signac.get_project().root_directory(),
                     'xyz_files',
                     job.doc.atoms_filename)))
    # We're assuming this needs to be run in the workspace directory
    # because of the sparseX files that go along with the potential
    args.append('param_filename=potential.xml')
    args.append('E')
    if do_forces:
        args.append('F')
    args.append('timing')
    return [cmd, ] + args


@FlowProject.label
def gap_fit_success(job):
    return job.isfile('potential.xml')


@FlowProject.operation
@FlowProject.pre(atoms_file_available)
@FlowProject.post(gap_fit_success)
@directives(omp_num_threads=18)
def fit_gap(job):
    cmdline = build_gap_fit_command_line(job)
    outfile_name = 'fit_gap.out'
    print(cmdline)
    with open(job.fn(outfile_name), 'a') as outfile:
        subprocess.run(cmdline, stdout=outfile)


@FlowProject.label
def has_energy_timings(job):
    return ('energy_ip_total_time_peratom_mean' in job.doc)

@FlowProject.label
def has_force_timings(job):
    return ('force_ip_total_time_peratom_mean' in job.doc)


def has_libatoms_exit_message(filename):
    """Check if the (potentially very large) file ends in the expected string"""
    fileend = "libAtoms::Finalise: Bye-Bye!\n"
    # Warning: Potential race condition(s)
    # (also implicitly assuming UTF-8 encoding)
    fsize = os.stat(filename).st_size
    with open(filename, 'r') as outfile:
        outfile.seek(fsize - len(fileend))
        return (outfile.read(len(fileend)) == fileend)


def has_energy_output(job):
    outfile_name = 'energies_output.out'
    if not job.isfile(outfile_name):
        return False
    else:
        return has_libatoms_exit_message(job, job.fn(outfile_name))


def has_force_output(job):
    outfile_name = 'forces_output.out'
    if not job.isfile(outfile_name):
        return False
    else:
        return has_libatoms_exit_message(job, job.fn(outfile_name))


@FlowProject.operation
@FlowProject.pre(gap_fit_success)
@FlowProject.post(has_energy_output)
@directives(omp_num_threads=1)
def time_quip_energies(job):
    old_dir = os.getcwd()
    quip_cmdline = build_quip_command_line(job)
    print(quip_cmdline)
    try:
        os.chdir(job.workspace())
        with open('energies_output.out', 'w') as outfile:
            subprocess.run(quip_cmdline, stdout=outfile)
    finally:
        os.chdir(old_dir)


@FlowProject.operation
@FlowProject.pre(gap_fit_success)
@FlowProject.post(has_force_output)
@directives(omp_num_threads=1)
def time_quip_forces(job):
    """Run QUIP timings with forces (and energies too)"""
    old_dir = os.getcwd()
    quip_cmdline = build_quip_command_line(job, do_forces=True)
    print(quip_cmdline)
    try:
        os.chdir(job.workspace())
        with open('forces_output.out', 'w') as outfile:
            subprocess.run(quip_cmdline, stdout=outfile)
    finally:
        os.chdir(old_dir)


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
            if re.match('TIMER: calc_connect', line):
                connect_timer.append(float(line.split()[7]))
            if re.match('TIMER: soap_calc', line):
                soap_timer.append(float(line.split()[7]))
            if re.match('TIMER: IPModel_GAP_Calc_gp_predict', line):
                gap_timer.append(float(line.split()[7]))
            if re.match('TIMER: IP_Calc', line):
                ip_timer.append(float(line.split()[7]))
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
    if len(connect_timer) != 2*n_configs:
        raise ValueError(error_message.format(number=len(connect_timer),
                                              time_type='calc_connect',
                                              filename=job.fn(filename)))
    for time_list, time_name in zip((soap_timer, gap_timer, ip_timer),
                                    ('soap_calc', 'gap_calc', 'IP_calc')):
        if len(time_list) != n_configs:
            raise ValueError(error_message.format(
                number=len(time_list),
                time_type=time_name,
                filename=job.fn(filename)
            )
    return {'natoms': np.array(natoms_list),
            'calc_connect': np.array(connect_timer),
            'soap': np.array(soap_timer),
            'gap': np.array(gap_timer),
            'ip_total': np.array(ip_timer)}


def store_timing_data(job, data, timing_key):
    """Store QUIP timing data and compute summary statistics"""
    for key, timer in data.items():
        datakey = '_'.join(timing_key, key)
        job.data[datakey] = timer
    natoms_list = data['natoms']
    del data['natoms']
    for key, timer in data.items():
        # Assuming calc_connect is called twice for each configuration
        if key == 'calc_connect':
            timer = np.sum(np.array(timer).reshape((-1, 2)), axis=1)
        job.doc['{:s}_{:s}_time_peratom_mean'.format(timing_key, key)] = np.mean(
            np.array(timer) / natoms_list)
        job.doc['{:s}_{:s}_time_peratom_min'.format(timing_key, key)] = np.min(
            np.array(timer) / natoms_list)
        job.doc['{:s}_{:s}_time_peratom_max'.format(timing_key, key)] = np.max(
            np.array(timer) / natoms_list)
        job.doc['{:s}_{:s}_time_peratom_std'.format(timing_key, key)] = np.std(
            np.array(timer) / natoms_list)


@FlowProject.operation
@FlowProject.pre(has_energy_output)
@FlowProject.post(has_energy_timings)
def process_energy_timings(job):
    timing_data = parse_quip_timings(job, 'energies_output.out')
    store_timing_data(job, timing_data, 'energy')


@FlowProject.operation
@FlowProject.pre(has_force_output)
@FlowProject.post(has_force_timings)
def process_force_timings(job):
    timing_data = parse_quip_timings(job, 'forces_output.out')
    store_timing_data(job, timing_data, 'force')


if __name__ == "__main__":
    FlowProject().main()
