#!/usr/bin/env python3
""" Task manager for running AlphaFold on Summit.

    Ingests a file of proteins and distributes work to process them among
    allocated dask workers.  The dask workers will be assigned to a single
    GPU and so any CPUs.

usage: alphafoldtaskmgr.py [-h] [--scheduler-timeout SCHEDULER_TIMEOUT] --scheduler-file SCHEDULER_FILE --input-file INPUT_FILE

AlphaFold task manager

optional arguments:
  -h, --help            show this help message and exit
  --scheduler-timeout SCHEDULER_TIMEOUT, -t SCHEDULER_TIMEOUT
                        dask scheduler timeout
  --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
  --input-file INPUT_FILE, -i INPUT_FILE
                        file containing proteins to process
"""
import argparse
import logging
import csv
from pathlib import Path
from time import time, sleep

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

from distributed import Client, LocalCluster, as_completed


def get_num_workers(client):
    """
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())


def read_input_file(input_file):
    """ Read text file containing proteins to be processed

    :param input_file: path to input file of proteins
    :return: list of proteins to be processed
    """
    input_file_path = Path(input_file)
    if not input_file_path.exists():
        logging.critical(f'{input_file} does not exist')
        raise RuntimeError(f'{input_file} does not exist')

    with input_file_path.open('r') as file_input:
        # We call strip(x) to nuke trailing newline
        return list(x.strip() for x in file_input)


def run_alphafold(protein):
    """ This is the task sent to the workers to run alphafold
    :param protein: is the protein to be processed for this task
    :returns: start and stop times as well as protein processed
    """
    start_time = time()

    # TODO swap in subprocess call to alphafold here
    sleep(3)

    stop_time = time()

    return start_time, stop_time, protein


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AlphaFold task manager')

    parser.add_argument('--scheduler-timeout', '-t', default=5000, type=int,
                        help='dask scheduler timeout')
    parser.add_argument('--scheduler-file', '-s', required=True,
                        help='dask scheduler file')
    parser.add_argument('--input-file', '-i', required=True,
                        help='file containing proteins to process')

    args = parser.parse_args()

    logging.info(f'Scheduler file: {args.scheduler_file}')
    logging.info(f'Scheduler timeout: {args.scheduler_timeout}')
    logging.info(f'Input file: {args.input_file}')

    # Slurp in all the proteins we have to process
    proteins = read_input_file(args.input_file)

    logging.info(f'Read {len(proteins)} proteins to process.')

    with Client(LocalCluster(),  # scheduler_file=args.scheduler_file,
                timeout=args.scheduler_timeout,
                name='alphafoldtaskmgr') as client:

        logging.info(f'Starting with {get_num_workers(client)} dask workers.')

        task_futures = client.map(run_alphafold, proteins)

        ac = as_completed(task_futures)

        for i, finished_task in enumerate(ac):
            start_time, stop_time, protein = finished_task.result()

            logging.info(f'{protein} processed in '
                         f'{(stop_time - start_time) / 60} minutes.')
            logging.info(f'{len(proteins) - i} proteins left')

        logging.info(f'Finished with {get_num_workers(client)} dask workers '
                     f'still active.')

    logging.info('Done.')
