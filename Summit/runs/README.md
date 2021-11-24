# Run scripts and support files

This directory contains scripts for running AlphaFold on Summit.

## Test lists of proteins

* `casp14_fm.lst` -- 32 proteins used to test Summit scripts
* `test.lst` -- essentially just the first two lines of `casp14_fm.lst` for even
  smaller tests

## Environment files

* `test-env-file` contains environment variables to tweak memory settings; used
  in singularity invocations

## alphataskmgr support

* `many_nodes.lsf` -- LSF script for alphataskmgr on Summit
* `common.env` -- common environment variables for `many_nodes.lsf`
* `spawn_worker.sh` -- shell script invoked by `many_nodes.lsf` to assigned 
  `dask` workers specific GPUs

# Running alphataskmgr on Summit

## Pre-conditions

* That the Summit conda environment 
  `/gpfs/alpine/bip198/proj-shared/mcoletti/conda` exists
  * It better exist
  * It contains all the python dependencies needed by `alphataskmgr.py`
  * It's a clone from the OLCF supplied python `open-ce` conda 
* You've cloned this repo on Summit
  * Hopefully in the gpfs since that's high speed

## To run it

You will need to update some files in the clone of this repo on Summit 
before submitting a job.

1. Change `many_nodes.lsf` for your run in the Summit clone of this repo
   1. Change the BSUB headings to desired job characteristics
   2. Fix paths 
      1. `RUN_DIR` needs to point where you expect to run the script to 
         capture the output
      2. The path for `common.env` needs to point to where you cloned the 
         repo on gpfs
      3. The path to call to `spwan_worker.sh` needs synced to where you 
         cloned this repo
   3. Point to the protein list
      1. The `--input-file` in the `alphafoldtaskmgr` invocation needs to 
         point to the list of proteins that you wish to process
      2. E.g., `/gpfs/alpine/world-shared/bif135/desulfovibrio/de_hildenborough.lst`
2. Change `spawn_worker.sh`
   1. `common.env` path needs fixed here, too
   2. `RUN_DIR`, too
