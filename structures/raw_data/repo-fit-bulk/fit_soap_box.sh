#!/bin/bash

export PATH=~/QUIP/build/linux_x86_64_gfortran_openmp:$PATH

AT_FILE=mebox-minimal-pbe-b1b.xyz
GP_FILE=gp-mebox-pbe-ts-mbd.xml

time teach_sparse at_file=$AT_FILE \
    gap={soap atom_sigma=0.5 l_max=8 n_max=8 cutoff=6.0 cutoff_transition_width=1.0 delta=0.01 add_species n_species=2 species_z={{1 6}} n_sparse=2000 covariance_type=dot_product sparse_method=cur_points zeta=4.0} \
    default_sigma={0.0001 0.002 1.0 1.0} sparse_jitter=1e-10 virial_parameter_name=none gp_file=$GP_FILE
