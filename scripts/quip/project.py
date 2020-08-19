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
                            "Slurm substitutions (e.g. %j) are accepted.")
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
    args.append('gp_file={:s}'.format(job.fn('potential.xml')))
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
    print(cmdline)
    subprocess.run(cmdline)


if __name__ == "__main__":
    FlowProject().main()
