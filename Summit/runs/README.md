# Run scripts and support files

This directory and sub-directories contain scripts for running AlphaFold on
Summit.

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
