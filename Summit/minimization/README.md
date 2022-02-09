
## SETTING UP FOR OPENMM/AMBER MINIMIZATION ON SUMMIT

First, install `miniconda`:
```bash
module load cuda/10.2.89 gcc/8.1.1
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-ppc64le.sh
bash Miniconda3-latest-Linux-ppc64le.sh -b -p miniconda
# Initialize your ~/.bash_profile
miniconda/bin/conda init bash
source ~/.bashrc
```

The `ppc64le` packages have been uploaded to the [`omnia`](https://anaconda.org/omnia) and [`conda-forge`](https://anaconda.org/conda-forge) channels:
```bash
# Add conda-forge and omnia to your channel list
conda config --add channels omnia --add channels conda-forge
# Update to conda-forge versions of packages
conda update --yes --all
```
### Installing OpenMM and other packages

```bash
# Create a new environment named 'openmm'
conda create -n openmm python==3.7
# Activate it
conda activate openmm
# Install the 'openmm' 7.4.0 dev package for ppc64le into this environment
conda install --yes -c omnia-dev/label/cuda101 openmm
conda install pdbfixer dask
# start an interactive job on a single node of SUMMIT
bsub -W 2:00 -nnodes 1 -P bip198 -alloc_flags gpudefault -Is /bin/bash
# test installation of OpenMM
python -m openmm.testInstallation
```

## USAGE: 

```bash
bsub many_nodes.sh	# after editing the bsub lines
```
