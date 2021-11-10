#!/bin/bash
#
# Batch submission script for testing post-AF minimization run on Summit
#
# Support issue-14 branch on PSP git repository
#
#BSUB -P BIP198
#BSUB -W 00:30
#BSUB -q debug
#BSUB -nnodes 1
#BSUB -alloc_flags gpudefault
#BSUB -J af_min
#BSUB -o testing.%J.out
#BSUB -e testing.%J.err

RUN_DIR=/gpfs/alpine/proj-shared/bip198/minimize_af/WP_010940344.1/
SCHEDULER_FILE=${RUN_DIR}/scheduler_file.json

source /autofs/nccs-svm1_home1/davidsonrb/Scripts/debugging/af_minimization/common.env

mkdir -p $RUN_DIR
cd $RUN_DIR

# Copy over the hosts allocated for this job so that we can later verify
# that all the allocated nodes were busy with the correct worker allocation.
cat $LSB_DJOB_HOSTFILE | sort | uniq > $LSB_JOBID.hosts

# We need to figure out the number of nodes to later spawn the workers
NUM_NODES=$(cat $LSB_DJOB_HOSTFILE | sort | uniq | wc -l)
export NUM_NODES=$(expr $NUM_NODES - 1)

##
## Start dask scheduler on an arbitrary couple of CPUs (more than one CPU to
## handle overhead of managing all the dask workers.)
##

# The scheduler doesn't need GPUs. We give it two CPUs and 8 total threads to
# handle the overhead of managing so many workers.
jsrun --smpiargs="off" --gpu_per_rs 0 --nrs 1 --tasks_per_rs 1 --cpu_per_rs 2 --rs_per_host 1  --latency_priority cpu-cpu -b none dask-scheduler --interface ib0 --no-dashboard --no-show --scheduler-file $SCHEDULER_FILE &

# Give the scheduler a chance to spin up.
sleep 5

##
## Start the dask-workers, which will be paired up to an individual GPU.  This
## bash script will manage the dask workers and GPU allocation for each Summit node.
##

# Now launch ALL the dask workers simultaneously.  They won't come up at the
# same time, though.
jsrun -h $RUN_DIR --smpiargs="off" --nrs $NUM_NODES -e individual --stdio_stdout ${RUN_DIR}/worker_out.%j.%h.%p --stdio_stderr ${RUN_DIR}/worker_error.%j.%h.%p --tasks_per_rs 1 --cpu_per_rs 39 --gpu_per_rs 6 --rs_per_host 1 -b none --latency_priority gpu-cpu /autofs/nccs-svm1_home1/davidsonrb/Scripts/debugging/af_minimization/spawn_worker.sh &

echo Waiting for workers

# Hopefully long enough for some workers to spin up and wait for work
sleep 30

# Run the client task manager; like the scheduler, this just needs a single
# core to noodle away on.
jsrun  --smpiargs="off" -e individual --stdio_stdout ${RUN_DIR}/alphafold_stdout.%j.%h.%p --stdio_stderr ${RUN_DIR}/alphafold_stderr.%j.%h.%p --nrs 1 --cpu_per_rs 1 --gpu_per_rs 0 --rs_per_host 1 --latency_priority cpu-cpu -b none python3 -u /autofs/nccs-svm1_home1/davidsonrb/Scripts/debugging/af_minimization/minimization_taskmgr.py --scheduler-file $SCHEDULER_FILE --input-file /autofs/nccs-svm1_home1/davidsonrb/Scripts/debugging/af_minimization/input_file.lst

# We're done so kill the scheduler and worker processes
jskill all

echo Run finished.

