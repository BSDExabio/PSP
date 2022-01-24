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
import platform
import os
import stat

import pdbfixer
import openmm
import csv

from distributed import Client, as_completed, get_worker

# NOTE: hard coded variables for the moment.
RESTRAINT_SET = "non_hydrogen"  # or "c_alpha"
RELAX_EXCLUDE_RESIDUES = []     # fed to openmm minimizeEnergy; left empty by default


#######################################
### DASK RELATED FUNCTIONS
#######################################

def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s      %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)        
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


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
        return [x.rstrip() for x in file_input]


def will_restrain(atom: openmm.app.topology.Atom, rset: str) -> bool:
  """Returns True if the atom will be restrained by the given restraint set."""

  if rset == "non_hydrogen":
    return atom.element.name != "hydrogen"
  elif rset == "c_alpha":
    return atom.name == "CA"


def _add_restraints(
    system,             #system: openmm.System,
    reference_pdb,      #reference_pdb: openmm.PDBFile,
    stiffness,          #stiffness: openmm.unit.Unit,
    rset,               #rset: str,
    exclude_residues):  #exclude_residues: Sequence[int]):
  """Adds a harmonic potential that restrains the end-to-end distance."""
  
  assert rset in ["non_hydrogen", "c_alpha"]

  force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
  force.addGlobalParameter("k", stiffness)
  for p in ["x0", "y0", "z0"]:
    force.addPerParticleParameter(p)

  for i, atom in enumerate(reference_pdb.topology.atoms()):
    if atom.residue.index in exclude_residues:
      continue
    if will_restrain(atom, rset):
      force.addParticle(i, reference_pdb.positions[i])
  system.addForce(force)


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


def fix_protein(input_pdb_file, output_pdb_file = 'protonated.pdb', logger = None):
    """
    """
    logger.info(f'Checking model for any required fixes (missing hydrogens and other atoms, etc).')
    logger.info(f'        Loading {input_pdb_file} to check for missing atoms and add hydrogens.')
    fixer = pdbfixer.PDBFixer(pdbfile=open(input_pdb_file,'r'))
    
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    logger.info(f'        Saving {output_pdb_file}.')
    with open(output_pdb_file,'w') as save_file:
        openmm.app.pdbfile.PDBFile.writeFile(fixer.topology,fixer.positions,file=save_file)
    
    os.chmod(output_pdb_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)

    return output_pdb_file


def prep_protein(pdb_file, restraint_set = "", exclude_residues = [], forcefield = "amber99sb.xml", restraint_stiffness = 10.0, platform = 'CUDA', energy_units = openmm.unit.kilocalories_per_mole, length_units = openmm.unit.angstroms, logger = None):
    """
    """
    
    logger.info(f'Preparing the simulation engine:')
    logger.info(f'        Loading {pdb_file} to create the OpenMM simulation object.')
    # load pdb file into an openmm Topology and coordinates object
    pdb = openmm.app.pdbfile.PDBFile(pdb_file)

    # set the FF and constraints objects
    logger.info(f'        Using {forcefield}.')
    force_field = openmm.app.forcefield.ForceField(forcefield)
    
    # prepare the restraints/constraints for the system
    logger.info(f'        Building HBond constraints as well as restraints on {restraint_set}.')
    constraints = openmm.app.HBonds # NOTE: check that bond lengths are good to begin with...
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    # prep and add restraints to the simulation system
    stiffness = restraint_stiffness * energy_units / (length_units**2)
    if stiffness > 0. * energy_units / (length_units**2):
        # code up the _add_restraints function
        _add_restraints(system, pdb, stiffness, restraint_set, exclude_residues)
    
    logger.info(f'        Creating the OpenMM simulation object.')
    # required to set this for prepping the simulation object
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)    # hard set because we won't be using it; still necessary to define for the creation of the simulation object
    # determine what hardware will be used to perform the calculations
    platform = openmm.Platform.getPlatformByName(platform)
    # prep the simulation object
    simulation = openmm.app.Simulation(pdb.topology, system, integrator, platform)
    # set the atom positions for the simulation's system's topology
    simulation.context.setPositions(pdb.positions)
    
    return simulation


def run_minimization(simulation,out_file_name,max_iterations = 0,energy_tolerance = 2.39,fail_attempts=100,energy_units = openmm.unit.kilocalories_per_mole,length_units = openmm.unit.angstroms, logger = None):
    """
    """
    tolerance = energy_tolerance * energy_units
    attempts = 0
    minimized = False
    logger.info(f'Running the minimization:')

    # grab initial energies and positions
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    einit = state.getPotentialEnergy().value_in_unit(energy_units)
    posinit = state.getPositions(asNumpy=True).value_in_unit(length_units)
    logger.info(f'        Starting energy: {einit} kcal mol-1')
    
    # attempt to minimize the structure
    while not minimized and attempts < fail_attempts:
        ## running minimization
        #simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
        #
        ## return energies and positions
        #state = simulation.context.getState(getEnergy=True, getPositions=True)
        #efinal = state.getPotentialEnergy().value_in_unit(energy_units)
        #positions = state.getPositions(asNumpy=True).value_in_unit(length_units)
        #logger.info(f'        Final energy: {efinal} kcal mol-1')
        #
        ## saving the final structure to a pdb
        #openmm.app.pdbfile.PDBFile.writeFile(simulation.topology,positions,file=open(out_file_name + '_%02d.pdb'%(attempts-1),'w'))
        #minimized = True
        
        attempts += 1
        try:
            # running minimization
            simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
            
            # return energies and positions
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            efinal = state.getPotentialEnergy().value_in_unit(energy_units)
            positions = state.getPositions(asNumpy=True).value_in_unit(length_units)
            logger.info(f'        Final energy: {efinal} kcal mol-1')
             
            # saving the final structure to a pdb
            with open(out_file_name + '_min_%02d.pdb'%(attempts-1),'w') as out_file:
                openmm.app.pdbfile.PDBFile.writeFile(simulation.topology,positions,file=out_file)
                os.chmod(out_file_name + '_min_%02d.pdb'%(attempts-1), stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
            minimized = True
        except Exception as e:
            logger.info(e)
    
    # update for more energy logging
    # calculate rmsd between pos and posinit
    logger.info(f"        dE = {efinal - einit} kcal mol^{-1}")
    
    if not minimized:
        raise ValueError(f"Minimization failed after {fail_attempts} attempts.")
    
    return out_file_name + '_min_%02d.pdb'%(attempts-1)


def run_pipeline(pdb_file, restraint_set = 'non_hydrogen', relax_exclude_residues = []):
    """
    """
    full_start_time = time.time()
    path_breakdown = pdb_file.split('/')
    path = pdb_file.split(path_breakdown[-1])[0]   # grabbing working dir path by removing the file name
    model_descriptor = path_breakdown[-1][:-4]     # grabbing a good file naming descriptor 
    min_logger = setup_logger('minimization_logger', path + model_descriptor + '_min.log')  # setting up the individual run's logging file
    
    # load pdb file and add missing atoms (mainly hydrogens)
    try:
        start = time.time()
        pdb_file = fix_protein(pdb_file,output_pdb_file = path + model_descriptor + '_protonated.pdb', logger = min_logger)
        min_logger.info(f'Finished preparing the protein model for minimization; took {time.time() - start} secs.')
    except:
        min_logger.exception(f"Preparation of the protein failed. Killing this worker's task.")
        clean_logger(min_logger)
        os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        
        return full_start_time, time.time(), f'{pdb_file} failed on protein preparation.'

    # send protein into an OpenMM pipeline, preparing the protein for a simulation
    try:
        start = time.time()
        simulation = prep_protein(pdb_file,restraint_set, exclude_residues = relax_exclude_residues, logger = min_logger)
        min_logger.info(f'Finished preparing the simulation engine; took {time.time() - start} secs.')
    except:
        min_logger.exception(f"Preparation of the simulation engine failed. Killing this worker's task.")
        clean_logger(min_logger)
        os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        
        return full_start_time, time.time(), f'{pdb_file} failed on simulation preparation.'

    # run the minimization protocol and output minimized structure
    try:
        start = time.time()
        final_pdb = run_minimization(simulation, path + model_descriptor, logger = min_logger)
        min_logger.info(f'Finished running the minimization calculation; took {time.time() - start} secs.')
    except:
        min_logger.exception(f"Simulation failed. Killing this worker's task.")
        clean_logger(min_logger)
        os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        
        return full_start_time, time.time(), f'{pdb_file} failed on simulation.'
        
    clean_logger(min_logger)
    os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)

    worker = get_worker()
    return platform.node(), worker.id, full_start_time, time.time(), final_pdb
    

#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments
    parser = argparse.ArgumentParser(description='Post-AF energy minimization task manager')
    parser.add_argument('--scheduler-timeout', '-t', default=5000, type=int, help='dask scheduler timeout')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-file', '-i', required=True, help='file containing proteins to process')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    args = parser.parse_args()

    # setting up the main logger
    main_logger = setup_logger('tskmgr_logger','tskmgr.log')
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Scheduler timeout: {args.scheduler_timeout}')
    main_logger.info(f'Input file: {args.input_file}')

    # create list of strings pointing to pdb files to be energy minimized
    proteins = read_input_file(args.input_file)
    main_logger.info(f'Read {len(proteins)} proteins to process.')

    # setting up timing log file
    timings_file = open(args.timings_file, 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','protein'])
    timings_csv.writeheader()

    # starting dask client
    client = Client(scheduler_file=args.scheduler_file,timeout=args.scheduler_timeout,name='energymintaskmgr')
    main_logger.info(f'Client information: {client}')
    NUM_WORKERS = get_num_workers(client)
    main_logger.info(f'Starting with {NUM_WORKERS} dask workers.')

    # waiting for workers...
    wait_start = time.time()
    client.wait_for_workers(n_workers=NUM_WORKERS)
    main_logger.info(f'Waited for {NUM_WORKERS} workers took {time.time() - wait_start} sec')
    workers_info = client.scheduler_info()['workers']
    connected_workers = len(workers_info)
    main_logger.info(f'{connected_workers} workers connected')

    # do the thing.
    task_futures = client.map(run_pipeline,proteins, restraint_set = RESTRAINT_SET, relax_exclude_residues = RELAX_EXCLUDE_RESIDUES) 

    # gather results
    ac = as_completed(task_futures)
    for i, finished_task in enumerate(ac):
        hostname, worker_id, start_time, stop_time, protein = finished_task.result()
        if 'failed' in protein:
            main_logger.info(f'{protein}')
            main_logger.info(f'{len(proteins) - i - 1} proteins left')
            append_timings(timings_csv,hostname,worker_id,start_time,stop_time,protein)
        else:
            main_logger.info(f'{protein} processed in {(stop_time - start_time) / 60.} minutes.')
            main_logger.info(f'{len(proteins) - i - 1} proteins left')
            append_timings(timings_csv,hostname,worker_id,start_time,stop_time,protein)

    # closing log files and shutting down the cluster
    timings_file.close()
    os.chmod(args.timings_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)
    os.chmod('tskmgr.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
    workers_list = list(workers_info)
    disconnect(client,workers_list)

