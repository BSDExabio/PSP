Scripts associated with the task 
[`issue-10`](https://github.com/BSDExabio/PSP/issues/10):

> Mu 12:59 PM
use model_[1-5]_ptmis preferred. It has extra features that are useful for applications. Also, please try --preset=casp14 instead of --preset=reduced_dbs, this option will increase the run time dramatically, so test this second option later and it may improve the model quality a little.

> I may want to add --preset as a command line argument to support quickly 
> changing presets.

I also added support for `--feature-dir`.

Since these changes are just for adding support for two command line arguments,
will dial down this run to a single node for 30 minutes in the `debug` queue.
I.e., run it just long enough to ensure the command line arguments are being
properly read in.
