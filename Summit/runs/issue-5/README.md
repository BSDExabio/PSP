Scripts associated with the task 
[`issue-5`](https://github.com/BSDExabio/PSP/issues/5), which is essentially
setting up a `dask`-based distributed task management system for AlphaFold
on Summit.

In this case, we're going to use a single Summit node with six `dask` workers,
one for each GPU.  We'll be using the test protein files.
