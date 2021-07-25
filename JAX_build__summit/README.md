# JAX on Summit

The purpose of this is to track the WIP JAX build on Summit. i

JAX uses Bazel for its build system. This directory contains:

* How we built Bazel on Summit (Bazel is required to build JAX)
* Toolchain changes to Bazel (enables the ability to  use Summit compiler modules)
* Current issues with building JAX on Summit
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
startup --output_user_root=/gpfs/alpine/<your_path>  # anywhere on gpfs is good; default is in /tmp which eats cgroup mem limit
build --jobs 2 --local_ram_resources=HOST_RAM*0.04 # throttles Bazels resource usage. Has been seen to help with JAX build.
test --jobs 2
```
### What to change in JAX/Bazel setup

This is all a WIP as we figure out how to work with Bazel.

* Create alternative toolchain to enable Bazel's access to Summit compiler modules (initial idea outlined in the JAX issue above).

### Executing build

Will be adding more here soon.
