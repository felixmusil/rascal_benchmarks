#!/bin/bash

export PATH=~/QUIP/build/linux_x86_64_gfortran_openmp:$PATH

AT_FILE=mebox-minimal-pbe0-mbdint.xyz
GP_FILE=gp-mbd-pbe0-soap.xml

time teach_sparse at_file=$AT_FILE core_param_file=dispts_quip_params.xml core_ip_args={Potential xml_label=ts calc_args={hirshfeld_vol_name=hirshfeld_avg_volume}} e0=0.0 \
    gap={soap atom_sigma=0.5 l_max=8 n_max=8 cutoff=5.0 cutoff_transition_width=1.0 delta=0.001 add_species n_species=2 species_z={{1 6}} n_sparse=2000 covariance_type=dot_product sparse_method=cur_points zeta=4.0} \
    default_sigma={0.0001 1.0 1.0 1.0} sparse_jitter=1e-10 gp_file=$GP_FILE
