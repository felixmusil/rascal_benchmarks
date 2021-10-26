{% extends "slurm.sh" %}
{% block project_header %}

module load gcc/8.3.0
module load openblas/0.3.6-openmp

/home/veit/miniconda3/bin/activate

export PATH=$PATH:/home/veit/QUIP/build/linux_x86_64_gfortran_openmp/
# Prevent mis-detection as lfs (and subsequent failure)
export SIGNAC_FLOW_ENVIRONMENT=HelvetiosEnvironment

echo $HOSTNAME

{% endblock %}

