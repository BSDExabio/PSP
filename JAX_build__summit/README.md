# JAX on Summit

The purpose of this is to track the WIP JAX build on Summit.

JAX uses Bazel for its build system. This directory contains:

* How we built Bazel on Summit (Bazel is required to build JAX)
* Toolchain changes to Bazel (enables the ability to  use Summit compiler modules)
* List of issues we have opened with other projects

Open Jax issue: https://github.com/google/jax/issues/4493

## Installing Bazel on Summit

Run ./compiles.sh on a batch or compute node. The cgroup limits on login node tend to kill the build. Error may not indicate that clearly. 
```
module purge
module load gcc/4.8.5 go/1.11.5 python/3.7.0 #needs a python3 interpreter which isn't available in /usr/bin until the RHEL8 upgrade. I thought it needed go, but maybe not? Loading that anyway.
wget https://github.com/bazelbuild/bazel/releases/download/4.1.0/bazel-4.1.0-dist.zip
unzip bazel-4.1.0-dist.zip
./compile.sh
```
The below is useful knowledge. It has been seen to help with building Bazel itself, but also help with projects that need to be built by Bazel.

```
$ cat ~/.bazelrc 
build --jobs 2 --local_ram_resources=HOST_RAM*0.04 # throttles Bazels resource usage. Has been seen to help with JAX build.
test --jobs 2
```
### What to change in JAX/Bazel setup

This is all a WIP as we figure out how to work with Bazel.

* Create alternative toolchain to enable Bazel's access to Summit compiler modules (initial idea outlined in the JAX issue above).

### Setup

#### Using GCC 7.4.0

* Load modules
```
module purge
module load python/3.7.0-anaconda3-5.3.0 gcc/7.4.0 cuda/10.2.89
```
* Clone JAX and use 'main' branch (ppc changes)
```
[rprout@login2.summit jax]$ pwd
/gpfs/alpine/stf007/scratch/rprout/jax
[rprout@login2.summit jax]$ git branch -v
* main 1db53b1 Merge pull request #7365 from hawkinsp:ppc
```
* Create 'toolchain' directory at the top of the JAX directory
* Add these files to the toolchain directory:
  * cc_toolchain_config.bzl https://gist.github.com/proutrc/54435117921a1fd581c4986448767b84
  * BUILD https://gist.github.com/proutrc/ba20b7f5d2c4b6ae7e26cfe4afdea44e
```
[rprout@login2.summit jax]$ cd toolchain/
[rprout@login2.summit toolchain]$ pwd
/gpfs/alpine/stf007/scratch/rprout/jax/toolchain
[rprout@login2.summit toolchain]$ ls
BUILD  cc_toolchain_config.bzl
```
* Throttle Bazel with this entry in your ~/.bazelrc file. 
  * NOTE: Without this, the build will spit out `undeclared inclusion(s)` errors (pointing at the GCC .h files). These errors are a bit of a mystery at the moment,           but they appear often when trying to build JAX with Bazel. In theory, the `cc_toochain_config.bzl` should fix this. It does fix it mostly, but only when             we throttle Bazel like this. Might be more to play with here.
```
$ cat ~/.bazelrc 
build --jobs 2 --local_ram_resources=HOST_RAM*0.04 # throttles Bazels resource usage. Has been seen to help with JAX build.
test --jobs 2
```
* Establish Bazel disk cache (Pick a place you prefer). Add this to your `~/.bazelrc` file too.
```
build --disk_cache=/gpfs/alpine/stf007/scratch/rprout/bazel-cache/
```
NOTE: These `~/.bazelrc` options get inheritied at build time.
### Running build

* Setup TEST_TMPDIR environment variable. These is necessary in order to move Bazel's build I/O off NFS (NFS is problematic for Bazel on many levels it seems).
  * NOTE: GPFS is used here, but /tmp is another valid option.
```
export TEST_TMPDIR=/gpfs/alpine/stf007/scratch/rprout/bazel-tmp/
```
* It seems you can never run this command enough (when in doubt clean the Bazel cache, etc.). 
  * NOTE: You need to know where your Bazel binary is. 
```
/gpfs/alpine/stf007/scratch/rprout/bazel-4.1.0/output/bazel clean --expunge
rm -rf /gpfs/alpine/stf007/scratch/rprout/bazel-cache/
rm -rf /gpfs/alpine/stf007/scratch/rprout/bazel-tmp/
```
* Run the build (points at our toolchain files)
```
python build/build.py --bazel_path=/gpfs/alpine/stf007/scratch/rprout/bazel-4.1.0/output/bazel --noenable_mkl_dnn --enable_cuda --cuda_path $OLCF_CUDA_ROOT --cudnn_path $OLCF_CUDA_ROOT --target_cpu=ppc --bazel_options=--host_crosstool_top=//toolchain:ppc --bazel_options=--crosstool_top=//toolchain:ppc --bazel_options=--cpu=ppc
```
* Should get this reproducible error : https://gist.github.com/proutrc/a4213af07897d521512d9d652b5ed0ad 

#### Using GCC 8.1.1

* Using 8.1.1 seems to get us passed the above `error: 'vec_neg' was not declared in this scope` with 7.4.0.
* BUT, it starts giving us `undeclared inclusion(s)` again!! Seeming further in the build process.
* To use GCC 8.1.1 you will want to do these things:
  * Replace 7.4.0 module with 8.1.1
  * update `cc_toolchain_config.bzl` (both `tool_paths` and `inclusions`)
    * I have them commented out in the `cc_toolchain_config.bzl` link above.   
  * Clean the Bazel environment. Run the `--expunge` command and remove bazel-[cache/tmp], as noted above. 
