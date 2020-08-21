"""Run QUIP benchmarks"""

import os
import errno
import subprocess

import signac
import flow
from flow import FlowProject, directives


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
    return ('time_total_eonly' in job.doc)


@FlowProject.operation
@FlowProject.pre(gap_fit_success)
@FlowProject.post(has_energy_timings)
@directives(omp_num_threads=1)
def time_quip_energies(job):
    old_dir = os.getcwd()
    quip_cmdline = build_quip_command_line(job)
    try:
        os.chdir(job.workspace())
        with open('energies_output.out', 'a') as outfile:
            subprocess.run(quip_cmdline, stdout=outfile)
    finally:
        os.chdir(old_dir)


if __name__ == "__main__":
    FlowProject().main()
