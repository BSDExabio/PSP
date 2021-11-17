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

import pdbfixer
from simtk import openmm
from simtk import unit
from simtk.openmm import app as openmm_app
from simtk.openmm.app.internal.pdbstructure import PdbStructure

#from rich.console import Console
#from rich.table import Table
#from rich import print
#from rich import pretty
#pretty.install()
#from rich.traceback import install
#install()
#from rich.logging import RichHandler
#rich_handler = RichHandler(rich_tracebacks=True,
#                           markup=True)
#logging.basicConfig(level='INFO', format='%(message)s',
#                    datefmt="[%Y/%m/%d %H:%M:%S]",
#                    handlers=[rich_handler])

from distributed import Client, as_completed

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


def will_restrain(atom: openmm_app.Atom, rset: str) -> bool:
  """Returns True if the atom will be restrained by the given restraint set."""

  if rset == "non_hydrogen":
    return atom.element.name != "hydrogen"
  elif rset == "c_alpha":
    return atom.name == "CA"


def _add_restraints(
    system,             #system: openmm.System,
    reference_pdb,      #reference_pdb: openmm_app.PDBFile,
    stiffness,          #stiffness: unit.Unit,
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
    openmm_app.pdbfile.PDBFile.writeFile(fixer.topology,fixer.positions,file=open(output_pdb_file,'w'))
    return output_pdb_file


def prep_protein(pdb_file, restraint_set = "", exclude_residues = [], forcefield = "amber99sb.xml", restraint_stiffness = 10.0, platform = 'CUDA', energy_units = unit.kilocalories_per_mole, length_units = unit.angstroms, logger = None):
    """
    """
    
    logger.info(f'Preparing the simulation engine:')
    logger.info(f'        Loading {pdb_file} to create the OpenMM simulation object.')
    # load pdb file into an openmm Topology and coordinates object
    pdb = openmm_app.PDBFile(pdb_file)

    # set the FF and constraints objects
    logger.info(f'        Using {forcefield}.')
    force_field = openmm_app.ForceField(forcefield)
    
    # prepare the restraints/constraints for the system
    logger.info(f'        Building HBond constraints as well as restraints on {restraint_set}.')
    constraints = openmm_app.HBonds # NOTE: check that bond lengths are good to begin with...
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
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    # set the atom positions for the simulation's system's topology
    simulation.context.setPositions(pdb.positions)
    
    return simulation


def run_minimization(simulation,out_file_name,max_iterations = 0,energy_tolerance = 2.39,fail_attempts=100,energy_units = unit.kilocalories_per_mole,length_units = unit.angstroms, logger = None):
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
        # running minimization
        simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
        
        # return energies and positions
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        efinal = state.getPotentialEnergy().value_in_unit(ENERGY)
        positions = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
        logger.info(f'        Final energy: {efinal} kcal mol-1')
        
        # saving the final structure to a pdb
        openmm_app.pdbfile.PDBFile.writeFile(simulation.topology,positions,file=open(out_file_name + '_%02d.pdb'%(attempts-1),'w'))
        minimized = True
        
        #attempts += 1
        #try:
        #    # running minimization
        #    simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
        #    
        #    # return energies and positions
        #    state = simulation.context.getState(getEnergy=True, getPositions=True)
        #    efinal = state.getPotentialEnergy().value_in_unit(ENERGY)
        #    positions = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
        #    
        #    # 
        #    openmm_app.pdbfile.PDBFile.writeFile(simulation.topology,positions,file=open(out_file_name + '_%02d.pdb'%(attempts-1),'w'))
        #    minimized = True
        #except Exception as e:
        #    logging.info(e)
    
    # update for more energy logging
    # calculate rmsd between pos and posinit
    logger.info(f"        dE = {efinal - einit} kcal mol^{-1}")
    
    if not minimized:
        raise ValueError(f"Minimization failed after {fail_attempts} attempts.")


def run_pipeline(pdb_file, restraint_set = 'non_hydrogen', relax_exclude_residues = []):
    """
    """
    path_breakdown = pdb_file.split('/')
    path = pdb_file.split(path_breakdown[-1])[0]   # grabbing working dir path by removing the file name
    model_descriptor = path_breakdown[-1].split('unrelaxed_')[-1][:-4]   # grabbing a good file naming descriptor 
    min_logger = setup_logger('minimization_logger', path + model_descriptor + '.log')  # setting up the individual run's logging file
    
    # load pdb file and add missing atoms (mainly hydrogens)
    start = time.time()
    pdb_file = fix_protein(pdb_file,output_pdb_file = path + model_descriptor + '_protonated.pdb', logger = min_logger)
    min_logger.info(f'Finished preparing the protein model for minimization; took {time.time() - start} secs.')

    # send protein into an OpenMM pipeline, preparing the protein for a simulation
    start = time.time()
    simulation = prep_protein(pdb_file,restraint_set, exclude_residues = relax_exclude_residues, logger = min_logger)
    min_logger.info(f'Finished preparing the simulation engine; took {time.time() - start} secs.')

    # run the minimization protocol and output minimized structure
    start = time.time()
    run_minimization(simulation, path + model_descriptor + '_min.pdb', logger = min_logger)
    min_logger.info(f'Finished running the minimization calculation; took {time.time() - start} secs.')


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

    # setting up the main logger
    main_logger = setup_logger('tskmgr_logger','tskmgr.log')
    main_logger.info('Starting dask pipeline and setting up logging.')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Scheduler timeout: {args.scheduler_timeout}')
    main_logger.info(f'Input file: {args.input_file}')

    # create list of strings pointing to pdb files to be energy minimized
    proteins = read_input_file(args.input_file)
    main_logger.info(f'Read {len(proteins)} proteins to process.')

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
    #task_futures = client.map(run_hello_world,list(range(NUM_WORKERS)))
    #task_futures = client.map(run_minimization,proteins)
    task_futures = client.map(run_pipeline,proteins, restraint_set = RESTRAINT_SET, relax_exclude_residues = RELAX_EXCLUDE_RESIDUES) 

    # gather results
    ac = as_completed(task_futures)
    for i, finished_task in enumerate(ac):
        start_time, stop_time, protein = finished_task.result()
        with open('out_%d.txt'%(i),'a') as out:
            out.write('total time spent on this: %f'%(stop_time-start_time))

        main_logger.info(f'{protein} processed in {(stop_time - start_time) / 60} minutes.')
        main_logger.info(f'{len(proteins) - i - 1} proteins left')


    # shutting down the cluster
    main_logger.info(f'Shutting down the cluster')
    workers_list = list(workers_info)
    disconnect(client,workers_list)
    main_logger.info('Done.')


    #with Client(scheduler_file=args.scheduler_file,timeout=args.scheduler_timeout,name='energymintaskmgr') as client:
    #    main_logger.info(f'Starting with {get_num_workers(client)} dask workers.')
    #    # setting the client task
    #    try:
    #        task_futures = client.map(run_minimization, proteins)
    #    except:
    #        print('something went wrong with the task_futures lines')
    #    
    #    # pushing tasks and watching them complete
    #    ac = as_completed(task_futures)
    #    for i, finished_task in enumerate(ac):
    #        start_time, stop_time, protein = finished_task.result()
    #        main_logger.info(f'{protein} processed in {(stop_time - start_time) / 60} minutes.')
    #        main_logger.info(f'{len(proteins) - i - 1} proteins left')
    #    
    #    main_logger.info(f'Finished with {get_num_workers(client)} dask workers still active.')

    #main_logger.info('Done.')


