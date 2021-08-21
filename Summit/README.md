# PSP on Summit

THIS IS A WIP.

This folder provides a way to get started on Summit. We currently only run the DL portion of AlphaFold on Summit. Preprocesing will likely take place on Andes.

## The container

We built the container on an ORNL Raptor machine, which is also Power9. The container build mechanism was Podman. We then converted it to a Singularity runtime file. 

The container file we run is: alphafold1103.sif (not exactly sure where to put the container yet, other than the Summit filesystem for sharing, it is too large to load into Github like a normal file).

It is built for CUDA11.

## The setup

- alphafold: Contains alphafold codebase, in addtion to our own scripts.  
- alphafold_databases: Ccontains parameters for the run. Also contains database tar. 
- casp14: Contains casp14 list file and reduced database for testing.
- desulfovibrio: A larger set of 559 sequences in the listi file.
- run: Contains run scripts. 


## How to do a test run within an interactive job

- Clone this directory.
- Make necessary changes to run/run_af_stage2a_batch_summit.sh (change directory base to represent your setup).
- example jsrun command:`jsrun -n1 -g6 singularity exec --bind /gpfs/alpine:/gpfs/alpine --nv /gpfs/alpine/stf007/scratch/rprout/AlphaFold/alphafold1103.sif ./run_af_stage2a_batch_summit.sh /gpfs/alpine/stf007/scratch/rprout/AlphaFold/runs/casp14_fm.lst /gpfs/alpine/stf007/scratch/rprout/AlphaFold/test-output` 

NOTE: You need to provide target list file and output directory to the run script

