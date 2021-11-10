#!/usr/bin/env bash
#
# Start up the dask workers with associated GPUs for a single
# Summit node.  This is called from an LSF batch script and is for a single
# resource set.
#
# Should be executed:
# jsrun -n 1 -c 8 -g 1 spawn_worker.sh
#

# Get common setup for alphafoldtaskmgr runs
source /autofs/nccs-svm1_home1/davidsonrb/Scripts/debugging/af_minimization/common.env

# Where we will run this script and expect output
RUN_DIR=/gpfs/alpine/proj-shared/bip198/minimize_af/WP_010940344.1/

# dask file for scheduler and workers to find each other
SCHEDULER_FILE=${RUN_DIR}/scheduler_file.json

echo "spawn worker python:" $(which python3)

# Grab the device numbers for all the local GPUs
gpus=$(nvidia-smi --list-gpus | cut -c5)


for gpu in $gpus; do
    env OMP_NUM_THREADS=4 SINGULARITYENV_SDL_VIDEODRIVER=offscreen CUDA_VISIBLE_DEVICES=$gpu \
      CARLA_PORT=$port dask-worker --nthreads 1 --nprocs 1 --interface ib0 \
      --no-dashboard --no-nanny --reconnect --scheduler-file ${SCHEDULER_FILE} &
done

wait
