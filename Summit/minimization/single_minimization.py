#!/usr/bin/env python3
"""
    Creating the model minimization process used in AF to be used for post-modeling AF runs
    
"""

#######################################
### PREAMBLE
#######################################
import argparse
import logging
import time
import numpy as np
import pdbfixer
from simtk import openmm
from simtk import unit
from simtk.openmm import app as openmm_app
from simtk.openmm.app.internal.pdbstructure import PdbStructure

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
### FUNCTIONS
#######################################
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


def fix_protein(input_pdb_file,output_pdb_file = 'protonated.pdb'):
    """
    """
    fixer = pdbfixer.PDBFixer(pdbfile=open(input_pdb_file,'r'))
    
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    openmm_app.pdbfile.PDBFile.writeFile(fixer.topology,fixer.positions,file=open(output_pdb_file,'w'))
    return output_pdb_file


def prep_protein(pdb_file,restraint_set,exclude_residues,forcefield="amber99sb.xml",restraint_stiffness = 10.0,platform='CUDA',energy_units = unit.kilocalories_per_mole,length_units = unit.angstroms):
    """
    """
    
    # load pdb file into an openmm Topology and coordinates object
    pdb = openmm_app.PDBFile(pdb_file)
    # set the FF and constraints objects
    force_field = openmm_app.ForceField(forcefield)
    constraints = openmm_app.HBonds # NOTE: check that bond lengths are good to begin with...
    # create the simulation ready system
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    
    # prep and add restraints to the simulation system
    stiffness = restraint_stiffness * energy_units / (length_units**2)
    if stiffness > 0. * energy_units / (length_units**2):
        # code up the _add_restraints function
        _add_restraints(system, pdb, stiffness, restraint_set, exclude_residues)
    
    # required to set this for prepping the simulation object
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)    # hard set because we won't be using it; still necessary to define for the creation of the simulation object
    # determine what hardware will be used to perform the calculations
    platform = openmm.Platform.getPlatformByName(platform)
    # prep the simulation object
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    
    positions = pdb.positions
    
    # set the atom positions for the simulation's system's topology
    simulation.context.setPositions(positions)
    
    return simulation


def run_minimization(simulation,out_file_name,max_iterations = 0,energy_tolerance = 2.39,fail_attempts=100,energy_units = unit.kilocalories_per_mole,length_units = unit.angstroms):
    """
    """

    
    tolerance = energy_tolerance * energy_units
    attempts = 0
    minimized = False

    # grab initial energies and positions
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    einit = state.getPotentialEnergy().value_in_unit(energy_units)
    posinit = state.getPositions(asNumpy=True).value_in_unit(length_units)
    
    # attempt to minimize the structure
    while not minimized and attempts < fail_attempts:
        attempts += 1
        try:
            # running minimization
            simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
            
            # return energies and positions
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            efinal = state.getPotentialEnergy().value_in_unit(ENERGY)
            positions = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
            
            # 
            openmm_app.pdbfile.PDBFile.writeFile(simulation.topology,positions,file=open(out_file_name[:-4] + '_%02d.pdb'%(attempts-1),'w'))
            minimized = True
        except Exception as e:
            print(e)
    
    # update for more energy logging
    # calculate rmsd between pos and posinit
    delta_e = efinal - einit 
    
    #print(f"{time_end} secs spent on minimization {attempts}\n{delta_e} kcal mol^{-1} energy change ({einit},{efinal})\n")
    
    if not minimized:
        raise ValueError(f"Minimization failed after {fail_attempts} attempts.")
    

def run_pipeline(pdb_file,restraint_set,relax_exclude_residues):
    """
    """

    # load pdb file and add missing atoms (mainly hydrogens)
    start = time.time()
    pdb_file = fix_protein(pdb_file)
    time_end = time.time() - start
    print(time_end)
    
    # send protein into an OpenMM pipeline, preparing the protein for a simulation
    start = time.time()
    simulation = prep_protein(pdb_file,restraint_set,relax_exclude_residues)
    time_end = time.time() - start
    print(time_end)
    
    # run the minimization protocol and output minimized structure
    start = time.time()
    run_minimization(simulation,pdb_file[:-4] + '_min.pdb')
    time_end = time.time() - start
    print(time_end)


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
        task_futures = client.map(run_pipeline, proteins, RESTRAINT_SET,RELAX_EXCLUDE_RESIDUES)
        
        # pushing tasks and watching them complete
        ac = as_completed(task_futures)
        for i, finished_task in enumerate(ac):
            start_time, stop_time, protein = finished_task.result()
            logging.info(f'{protein} processed in {(stop_time - start_time) / 60} minutes.')
            logging.info(f'{len(proteins) - i - 1} proteins left')
        
        logging.info(f'Finished with {get_num_workers(client)} dask workers still active.')

    logging.info('Done.')


