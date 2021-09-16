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
import subprocess
from pathlib import Path
from time import time

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

from distributed import Client, as_completed, get_worker


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
        # Each line contains the protein and the number of segments; we only
        # want the protein, so we grab that.
        return list(x.split()[0] for x in file_input)


def run_alphafold(protein, preset, feature_dir):
    """ This is the task sent to the workers to run alphafold
    :param protein: is the protein to be processed for this task
    :param preset: is the AlphaFold preset, and can be casp14 or reduced_dbs
    :param feature_dir: is the directory where AlphaFold features are found
    :returns: start and stop times as well as protein processed
    """
    import sys
    import platform

    start_time = time()

    args = []
    args.append('singularity')
    args.append('exec')
    args.append('--bind')
    args.append('/gpfs/alpine:/gpfs/alpine')
    args.append('--nv')
    args.append('--env-file')
    args.append(
        '/gpfs/alpine/bip198/proj-shared/mcoletti/PSP/Summit/runs/test-env-file')
    args.append('/gpfs/alpine/stf007/world-shared/subil/alphafold1103.sif')
    args.append('python3')
    args.append(
        '/gpfs/alpine/bip198/proj-shared/mcoletti/PSP/Summit/alphafold/run_alphafold_stage2a.py')
    args.append(f'--fasta_paths={protein}.fas')
    args.append(f'--preset={preset}')
    args.append(
        '--data_dir=/gpfs/alpine/world-shared/bif135/alphafold_onsummit/alphafold_databases/')
    args.append('--output_dir=.')

    # directory for protein feature description information files
    args.append(f'--feature_dir={feature_dir}')

    args.append(
        f'--model_names=model_1_ptm,model_2_ptm,model_3_ptm,model_4_ptm,model_5_ptm')
    # args.append('--benchmark') # this turns on computationally expensive benchmarking

    # TODO remove this reality check once we have a stable implementation
    sys.stdout.write('args: ' + str(args) + '\n')

    completed_process = subprocess.run(args=args, text=True,
                                       capture_output=True)

    sys.stdout.write(completed_process.stdout)
    sys.stdout.flush()
    sys.stderr.write(completed_process.stderr)
    sys.stderr.flush()

    if completed_process.returncode != 0:
        print(f'Process returned with code'
              f' {completed_process.returncode} for protein {protein}')

    stop_time = time()

    worker = get_worker()

    return platform.node(), worker.id, start_time, stop_time, protein


def append_timings(csv_writer, hostname, worker_id, start_time, stop_time,
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AlphaFold task manager')

    parser.add_argument('--scheduler-timeout', default=5000, type=int,
                        help='dask scheduler timeout')
    parser.add_argument('--preset', '-p',
                        default='casp14',
                        help='value for AlphaFold preset')
    parser.add_argument('--feature-dir',
                        default='/gpfs/alpine/world-shared/bif135/desulfovibrio/afold_fea',
                        help='Directory where protein features are found')
    parser.add_argument('--timings-file', '-t',
                        help='CSV file for protein processing timings')
    parser.add_argument('--scheduler-file', '-s', required=True,
                        help='dask scheduler file')
    parser.add_argument('--input-file', '-i', required=True,
                        help='file containing proteins to process')

    args = parser.parse_args()

    logging.info(f'Scheduler file: {args.scheduler_file}')
    logging.info(f'Scheduler timeout: {args.scheduler_timeout}')
    logging.info(f'Input file: {args.input_file}')
    logging.info(f'Timings file: {args.timings_file!s}')
    logging.info(f'Feature directory: {args.feature_dir}')
    logging.info(f'Preset: {args.preset}')

    # Slurp in all the proteins we have to process
    proteins = read_input_file(args.input_file)

    logging.info(f'Read {len(proteins)} proteins to process.')

    if args.timings_file:
        timings_file = open(args.timings_file, 'w')
        timings_csv = csv.DictWriter(timings_file,
                                     ['hostname',
                                      'worker_id',
                                      'start_time',
                                      'stop_time',
                                      'protein'])
        timings_csv.writeheader()

    start_time = time()

    with Client(scheduler_file=args.scheduler_file,
                timeout=args.scheduler_timeout,
                name='alphafoldtaskmgr') as client:

        logging.info(f'Starting with {get_num_workers(client)} dask workers.')

        task_futures = client.map(run_alphafold, proteins,
                                  preset=args.preset,
                                  feature_dir=args.feature_dir)

        ac = as_completed(task_futures)

        for i, finished_task in enumerate(ac):
            hostname, worker_id, start_time, stop_time, protein = finished_task.result()

            logging.info(f'{protein} processed in '
                         f'{(stop_time - start_time) / 60} minutes.')
            logging.info(f'{len(proteins) - i - 1} proteins left')

            if args.timings_file:
                append_timings(timings_csv,
                               hostname,
                               worker_id,
                               start_time,
                               stop_time,
                               protein)

        logging.info(f'Finished with {get_num_workers(client)} dask workers '
                     f'still active.')

    if args.timings_file:
        timings_file.close()

    logging.info(f'Finished in {(time() - start_time) / 60 / 60:0.2f} hours')
