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

## 

## `issue-N` directories

These directories correspond to support scripts for corresponding issues on
our [kanban board](https://github.com/BSDExabio/PSP/projects/1)
