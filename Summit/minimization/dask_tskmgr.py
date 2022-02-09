#!/usr/bin/env python3
""" Task manager for running the dask pipeline for post-AF minimization on Summit.

    Ingests a file of .pdb files and distributes work to process them among allocated dask workers. The dask workers will be assigned to a single GPU and associated CPUs.

    USAGE:
        python3 minimization_taskmgr.py [-h] [--scheduler-timeout SCHEDULER_TIMEOUT] --scheduler-file SCHEDULER_FILE --input-file INPUT_FILE --timings-file TIMINGS_FILE.csv

    INPUT:
        -h, --help      show this help message and exit
        --scheduler-timeout SCHEDULER_TIMEOUT, -t SCHEDULER_TIMEOUT
                        dask scheduler timeout; default: 5000 seconds
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
        --input-file INPUT_FILE, -i INPUT_FILE
                        file containing paths to pdb files to process
        --timings-file TIMINGS_FILE.csv, -ts TIMINGS_FILE.csv
                        CSV file for protein processing timings
        --working-dir /path/to/dir/, -wd /path/to/dir/
                        full path to the directory within which files will be written
        --script-path /path/to/dir/script.py, -sp /path/to/dir/script.py
                        full path to the python script to be run within the subprocess

    HARD CODED VARIABLES:
        RESTRAINT_SET: set to "non_hydrogen"; "c_alpha" is also acceptable
        RELAX_EXCLUDE_RESIDUES: set to an empty list; used to ignore residues when setting up restraints for atom positions.
        DIRECTORY: name of the directory within which output for each task will be written
"""

import time
import argparse
import platform
import os
import stat
import traceback
import numpy as np

import csv

import subprocess
from subprocess import CalledProcessError

import dask.config
from distributed import Client, Worker, as_completed, get_worker

import logging_functions

# NOTE: hard coded variables for the moment.
RESTRAINT_SET = "non_hydrogen"  # or "c_alpha"
RELAX_EXCLUDE_RESIDUES = []     # fed to openmm minimizeEnergy; left empty by default
DIRECTORY = 'relaxation'

#######################################
### DASK RELATED FUNCTIONS
#######################################

def get_num_workers(client):
    """ Get the number of active workers
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())


def disconnect(client, workers_list):
    """ Shutdown the active workers in workers_list
    :param client: active dask client
    :param workers_list: list of dask workers
    """
    client.retire_workers(workers_list, close_workers=True)
    client.shutdown()


def read_input_file(input_file):
    """ Read text file containing proteins to be processed
    :param input_file: str; a path to input file of proteins
    :return: list of strings that each point to a protein model to be processed
    """
    with open(input_file,'r') as file_input:
        # Each line contains the path to a AF model .pdb file
        pdb_list = [x.rstrip() for x in file_input]
        file_sizes = np.array([os.path.getsize(pdb) for pdb in pdb_list])
        return [pdb_list[i] for i in np.argsort(file_sizes)[::-1]]


def append_timings(csv_writer, file_object, hostname, worker_id, start_time, stop_time,
                   protein):
    """ append the protein timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param protein: that was processed
    """
    csv_writer.writerow({'hostname'  : hostname,
                         'worker_id' : worker_id,
                         'start_time': start_time,
                         'stop_time' : stop_time,
                         'protein'   : protein})
    file_object.flush()


def submit_pipeline(pdb_file,script,working_directory,restraint_set="non_hydrogen",relax_exclude_residues=[],save_directory=''):
    """
    """
    worker = get_worker()
    start_time = time.time()
    worker.logger.info(f'Starting to process {pdb_file} at {start_time}.')

    # FIXME temporarily commented out to fix logger error

    # try:
    #     completed_process = subprocess.run(['python3',script,pdb_file,restraint_set,save_directory],shell=False,capture_output=True,check=True,cwd=working_directory)
    #     final_pdb = completed_process.stdout.decode('utf-8').strip()
    #     worker.logger.info(f'Finished minimization of {pdb_file}; saved to {final_pdb}.')
    #     return platform.node(), worker.id, start_time, time.time(), final_pdb
    #
    # except CalledProcessError as e:
    #     worker.logger.error(f'Exception occurred. Return code {e.returncode}')
    #     worker.logger.error(f'Exception cmd: {e.cmd}')
    #     return platform.node(), worker.id, start_time, time.time(), f'failed to minimize {pdb_file}'
    #
    # except Exception as e:
    #     worker.logger.error(f'Odd exception {e} raised')
    #     worker.logger.error(f'Exception cmd: {e.cmd}')
    #     return platform.node(), worker.id, start_time, time.time(), f'failed to minimize {pdb_file}'


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Post-AF energy minimization task manager')
    parser.add_argument('--scheduler-timeout', '-t', default=5000, type=int, help='dask scheduler timeout')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-file', '-i', required=True, help='file containing proteins to process')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--working-dir', '-wd', required=True, help='path that points to the working directory for the output files')
    parser.add_argument('--script-path', '-sp', required=True, help='path that points to the script for the subprocess call')
    args = parser.parse_args()

    if args.working_dir[-1] != os.path.sep:
        args.working_dir += os.path.sep

    # set up the main logger file and list all relevant parameters.
    main_logger = logging_functions.setup_logger('tskmgr_logger',f'{args.working_dir}tskmgr.log')
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Scheduler timeout: {args.scheduler_timeout}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Working directory: {args.working_dir}')
    main_logger.info(f'Path to subprocess script: {args.script_path}')
    main_logger.info(f'Structure input file: {args.input_file}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')

    # create list of strings pointing to pdb files to be energy minimized.
    proteins = read_input_file(args.input_file)
    nProteins = len(proteins)
    main_logger.info(f'Read and sorted {nProteins} proteins for relaxation.')

    # set up timing log file.
    timings_file = open(args.timings_file, 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','protein'])
    timings_csv.writeheader()

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=args.scheduler_timeout,name='energymintaskmgr')
    main_logger.info(f'Client information: {client}')
    NUM_WORKERS = get_num_workers(client)
    main_logger.info(f'Starting with {NUM_WORKERS} dask workers.')

    # wait for workers.
    wait_start = time.time()
    client.wait_for_workers(n_workers=NUM_WORKERS)
    main_logger.info(f'Waited for {NUM_WORKERS} workers took {time.time() - wait_start} sec')
    workers_info = client.scheduler_info()['workers']
    connected_workers = len(workers_info)
    main_logger.info(f'{connected_workers} workers connected. Log files created.')

    client.upload_file('logging_functions.py')

    client.register_worker_plugin(logging_functions.WorkerLoggerPlugin())

    main_logger.info(f'Installed worker plugin.')

    # do the thing.
    task_futures = client.map(submit_pipeline,proteins,args.script_path,args.working_dir, restraint_set = RESTRAINT_SET, relax_exclude_residues = RELAX_EXCLUDE_RESIDUES, save_directory = DIRECTORY, pure=False)

    # gather results.
    ac = as_completed(task_futures)
    for i, finished_task in enumerate(ac):
        results = finished_task.result()
        hostname, worker_id, start_time, stop_time, protein = finished_task.result()
        if 'failed' in protein:
            main_logger.info(f'{protein}')
            main_logger.info(f'{nProteins - i - 1} proteins left')
            append_timings(timings_csv,timings_file,hostname,worker_id,start_time,stop_time,protein)
        else:
            main_logger.info(f'{protein} processed in {(stop_time - start_time) / 60.} minutes.')
            main_logger.info(f'{nProteins - i - 1} proteins left')
            append_timings(timings_csv,timings_file,hostname,worker_id,start_time,stop_time,protein)

    # close log files and shut down the cluster.
    timings_file.close()
    os.chmod(args.timings_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    logging_functions.clean_logger(main_logger)
    os.chmod('tskmgr.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
    client.shutdown()
