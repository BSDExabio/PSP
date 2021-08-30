#!/usr/bin/env python3
""" Task manager for running AlphaFold on Summit.

    Ingests a file of proteins and distributes work to process them among
    allocated dask workers.  The dask workers will be assigned to a single
    GPU and so any CPUs.
"""
import argparse
import logging
import csv

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
logging.basicConfig(level='DEBUG', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])

from distributed import Client, LocalCluster

def get_num_workers(client):
    """
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AlphaFold task manager')

    parser.add_argument('--scheduler-file', '-s', help='dask scheduler file')
    parser.add_argument('--scheduler-timeout', '-t', default=5000, type=int,
                        help='dask scheduler timeout')
    parser.add_argument('--input-file', '-i',
                        help='file containing proteins to process')

    args = parser.parse_args()

    logging.info(f'Scheduler file: {args.scheduler_file}')
    logging.info(f'Scheduler timeout: {args.scheduler_timeout}')
    logging.info(f'Input file: {args.input_file}')

    with Client(scheduler_file=args.scheduler_file,
                timeout=args.scheduler_timeout,
                name='alphafoldtaskmgr') as client:

        logging.info(f'Starting with {get_num_workers(client)} dask workers.')

        logging.info(f'Finished with {get_num_workers(client)} dask workers still active.')

    logging.info('Done.')
