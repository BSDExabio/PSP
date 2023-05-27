# PSP on Summit

This repo contains code use to run proteome-scale structure prediction pipelines using AlphaFold (v.2.0.) on the OLCF Summit supercomputer at the Oak Ridge National Laboratory.

We currently only run the DL portion of AlphaFold on Summit. Preprocesing (MSA generation, etc.) takes place on the Andes data analysis cluster.

## The container

We built the container on an ORNL Raptor machine, which is also Power9. The container build mechanism was Podman. We then converted it to a Singularity runtime file. See the followning papers for details: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9835405, https://ieeexplore.ieee.org/abstract/document/9652872.

The container file we run is: `alphafold1103.sif` (not exactly sure where to put the container yet, other than the Summit filesystem for sharing, it is too large to load into Github like a normal file).

It is built for CUDA11.

## The setup

- `alphafold`: Contains alphafold codebase, in addition to our own scripts.  
- `alphafold_databases/params`: Contains parameters for the models. 
- `casp14`: Contains casp14 list file and reduced database for testing (on GPFS at `/gpfs/alpine/world-shared/bif135/alphafold_onsummit/alphafold_test/casp14/af_reduced_db/`).
- `desulfovibrio:` A larger set of 559 sequences in the list file (on GPFS at `/gpfs/alpine/world-shared/bif135/desulfovibrio/afold_fea/`).
- `run:` Contains run scripts. 


## How to do a test run within an interactive job

- Clone this directory.
- Make necessary changes to run/run_af_stage2a_batch_summit.sh (change directory base to represent your setup).
- example jsrun command:

```jsrun -n1 -g6 singularity exec --bind /gpfs/alpine:/gpfs/alpine --nv /gpfs/alpine/stf007/scratch/rprout/AlphaFold/alphafold1103.sif ./run_af_stage2a_batch_summit.sh /gpfs/alpine/stf007/scratch/rprout/AlphaFold/runs/casp14_fm.lst /gpfs/alpine/stf007/scratch/rprout/AlphaFold/test-output``` 

NOTE: You need to provide target list file and output directory to the run script

