#!/bin/bash

sbatch ./model/submit_0.sh
sbatch ./model/submit_1.sh
sbatch ./model/submit_4.sh

sbatch ./spherical_expansion/submit_0.sh
sbatch ./spherical_expansion/submit_2.sh
sbatch ./spherical_expansion/submit_3.sh
sbatch ./spherical_expansion/submit_4.sh

sbatch ./spherical_invariants/submit_0.sh
sbatch ./spherical_invariants/submit_2.sh
sbatch ./spherical_invariants/submit_3.sh
sbatch ./spherical_invariants/submit_4.sh

# sbatch ./radial_integral/submit_0.sh

