#!/usr/bin/env python3
""" Task manager for running post-AF minimization on Summit. 

    Ingests a file of .pdb files and distributes work to process them among allocated dask workers. The dask workers will be assigned to a single GPU and associated CPUs. 

    USAGE: 
        python3 minimization_taskmgr.py [-h] [--scheduler-timeout SCHEDULER_TIMEOUT] --scheduler-file SCHEDULER_FILE --input-file INPUT_FILE

    INPUT: 
        -h, --help      show this help message and exit
        --scheduler-timeout SCHEDULER_TIMEOUT, -t SCHEDULER_TIMEOUT 
                        dask scheduler timeout; default: 5000 seconds
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
        --input-file INPUT_FILE, -i INPUT_FILE
                        file containing paths to pdb files to process

    HARD CODED VARIABLES:
        RESTRAINT_SET: set to "non_hydrogen"; "c_alpha" is also acceptable
        RELAX_EXCLUDE_RESIDUES: set to an empty list; used to ignore residues when setting up restraints for atom positions. 

"""

import subprocess
import time
import argparse
import logging

from rich.console import Console
from rich.table import Table
from rich import print
from rich import pretty
pretty.install()
from rich.traceback import install
install()
from rich.logging import RichHandler
rich_handler = RichHandler(rich_tracebacks=True,
                           markup=True)
logging.basicConfig(level='INFO', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])

from distributed import Client, as_completed

# NOTE: hard coded variables for the moment.
RESTRAINT_SET = "non_hydrogen"  # or "c_alpha"
RELAX_EXCLUDE_RESIDUES = []     # fed to openmm minimizeEnergy; left empty by default


#######################################
### DASK RELATED FUNCTIONS
#######################################

def get_num_workers(client):
    """
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())

def read_input_file(input_file):
    """ Read text file containing proteins to be processed
    :param input_file: str; a path to input file of proteins
    :return: list of strings that each point to a protein model to be processed
    """
    with open(input_file,'r') as file_input:
        # Each line contains the path to a AF model .pdb file
        return list(x for x in file_input)


def run_minimization(model):
    """
    """
    import sys

    start_time = time.time()

    # set up the job submit line to be used across all workers
    args = []
    args.append('jsrun')
    args.append('--smpiargs="none"')
    args.append('-n1')
    args.append('-a1')
    args.append('-c1')
    args.append('-g1')
    args.append('python3')
    args.append('/ccs/home/davidsonrb/Scripts/debugging/af_minimization/single_minimization.py')
    args.append(model)
    args.append(RESTRAINT_SET)
    args.append(RELAX_EXCLUDE_RESIDUES)

    sys.stdout.write('args: ' + str(args))

    completed_process = subprocess.run(args=args,text=True,capture_output=True)
    
    sys.stdout.write(completed_process.stdout)
    sys.stdout.flush()
    sys.stderr.write(completed_process.stderr)
    sys.stderr.flush()

    if completed_process.returncode != 0:
        print(f'Process returned with code'
              f' {completed_process.returncode} for protein {protein}')

    stop_time = time()

    return start_time, stop_time, protein


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments
    parser = argparse.ArgumentParser(description='Post-AF energy minimization task manager')
    parser.add_argument('--scheduler-timeout', '-t', default=5000, type=int, help='dask scheduler timeout')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-file', '-i', required=True, help='file containing proteins to process')
    args = parser.parse_args()

    logging.info(f'Scheduler file: {args.scheduler_file}')
    logging.info(f'Scheduler timeout: {args.scheduler_timeout}')
    logging.info(f'Input file: {args.input_file}')

    # create list of strings pointing to pdb files to be energy minimized
    proteins = read_input_file(args.input_file)

    logging.info(f'Read {len(proteins)} proteins to process.')

    # starting dask client
    with Client(scheduler_file=args.scheduler_file,timeout=args.scheduler_timeout,name='energymintaskmgr') as client:
        logging.info(f'Starting with {get_num_workers(client)} dask workers.')
        # setting the client task
        task_futures = client.map(run_minimization, proteins, RESTRAINT_SET,RELAX_EXCLUDE_RESIDUES)
        
        # pushing tasks and watching them complete
        ac = as_completed(task_futures)
        for i, finished_task in enumerate(ac):
            start_time, stop_time, protein = finished_task.result()
            logging.info(f'{protein} processed in {(stop_time - start_time) / 60} minutes.')
            logging.info(f'{len(proteins) - i - 1} proteins left')
        
        logging.info(f'Finished with {get_num_workers(client)} dask workers still active.')

    logging.info('Done.')


