#!/usr/bin/env python3

import pdbfixer
import openmm
import logging

import time
import argparse
import logging
import platform
import os
import stat
import traceback


#######################################
### PIPELINE FUNCTIONS
#######################################

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


def fix_protein(input_pdb_file, output_pdb_file = 'protonated.pdb'):
    """
    """
    logger_string = f'\nChecking model for any required fixes (missing hydrogens and other atoms, etc).\n        Loading {input_pdb_file} to check for missing atoms and add hydrogens.\n'
    
    fixer = pdbfixer.PDBFixer(pdbfile=open(input_pdb_file,'r'))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    logger_string += f'        Saving {output_pdb_file}.\n'
    with open(output_pdb_file,'w') as save_file:
        openmm.app.pdbfile.PDBFile.writeFile(fixer.topology,fixer.positions,file=save_file)
    
    os.chmod(output_pdb_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)

    del fixer

    return output_pdb_file, logger_string


def prep_protein(pdb_file, restraint_set = "", exclude_residues = [], forcefield = "amber99sb.xml", restraint_stiffness = 10.0, platform = 'CUDA', energy_units = openmm.unit.kilocalories_per_mole, length_units = openmm.unit.angstroms):
    """
    """
   
    logger_string = f'\nPreparing the simulation engine:\n        Loading {pdb_file} to create the OpenMM simulation object.\n'
    # load pdb file into an openmm Topology and coordinates object
    pdb = openmm.app.pdbfile.PDBFile(pdb_file)

    # set the FF and constraints objects
    logger_string += f'        Using {forcefield}.\n'
    force_field = openmm.app.forcefield.ForceField(forcefield)
    
    # prepare the restraints/constraints for the system
    logger_string += f'        Building HBond constraints as well as restraints on {restraint_set}.\n'
    constraints = openmm.app.HBonds # NOTE: check that bond lengths are good to begin with...
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    # prep and add restraints to the simulation system
    stiffness = restraint_stiffness * energy_units / (length_units**2)
    if stiffness > 0. * energy_units / (length_units**2):
        # code up the _add_restraints function
        _add_restraints(system, pdb, stiffness, restraint_set, exclude_residues)
   
    logger_string += f'        Creating the OpenMM simulation object.\n'
    # required to set this for prepping the simulation object
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)    # hard set because we won't be using it; still necessary to define for the creation of the simulation object
    # determine what hardware will be used to perform the calculations
    platform = openmm.Platform.getPlatformByName(platform)
    # prep the simulation object
    simulation = openmm.app.Simulation(pdb.topology, system, integrator, platform)
    # set the atom positions for the simulation's system's topology
    simulation.context.setPositions(pdb.positions)
   
    del pdb
    del force_field
    del constraints
    del system
    del stiffness
    del integrator
    del platform

    return simulation, logger_string


def run_minimization(simulation,out_file_name,max_iterations = 0,energy_tolerance = 2.39,fail_attempts=100,energy_units = openmm.unit.kilocalories_per_mole,length_units = openmm.unit.angstroms):
    """
    """
    tolerance = energy_tolerance * energy_units
    attempts = 0
    minimized = False

    # grab initial energies and positions
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    einit = state.getPotentialEnergy().value_in_unit(energy_units)
    posinit = state.getPositions(asNumpy=True).value_in_unit(length_units)
    
    logger_string = f'\nRunning the minimization:\n        Starting energy: {einit} kcal mol-1\n'
    
    # attempt to minimize the structure
    while not minimized and attempts < fail_attempts:
        attempts += 1
        try:
            # running minimization
            simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
            
            # return energies and positions
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            efinal = state.getPotentialEnergy().value_in_unit(energy_units)
            positions = state.getPositions(asNumpy=True).value_in_unit(length_units)
            logger_string += f'        Final energy: {efinal} kcal mol-1\n'
             
            # saving the final structure to a pdb
            with open(out_file_name + '_min_%02d.pdb'%(attempts-1),'w') as out_file:
                openmm.app.pdbfile.PDBFile.writeFile(simulation.topology,positions,file=out_file)
                os.chmod(out_file_name + '_min_%02d.pdb'%(attempts-1), stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
            minimized = True
        except Exception as e:
            logger_string += f'        Attempt {attempt}: {e}\n'
            logger.info(e)

    logger_string += f'        dE = {efinal - einit} kcal mol^{-1}\n'
    
    if not minimized:
        logger_string += f"Minimization failed after {fail_attempts} attempts.\n"
    
    return out_file_name + '_min_%02d.pdb'%(attempts-1), logger_string


### convert to if __name__ == '__main__': chunk of code
def run_pipeline(pdb_file, restraint_set = 'non_hydrogen', relax_exclude_residues = []):
    """
    """

    # load pdb file and add missing atoms (mainly hydrogens)
    try:
        start = time.time()
        pdb_file, logger_string = fix_protein(pdb_file,output_pdb_file = path + model_descriptor + '_protonated.pdb')
        logger_string += f'Finished preparing the protein model for minimization; took {time.time() - start} secs.\n'
        min_logger.info(logger_string)
    except Exception as e:
        traceback.print_exc()
        min_logger.exception(f"Preparation of the protein failed. Killing this worker's task.")
        clean_logger(min_logger)
        os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        
        return platform.node(), worker.id, full_start_time, time.time(), f'{pdb_file} failed on protein preparation.'

    # send protein into an OpenMM pipeline, preparing the protein for a simulation
    try:
        start = time.time()
        simulation, logger_string = prep_protein(pdb_file,restraint_set, exclude_residues = relax_exclude_residues)
        logger_string += f'Finished preparing the simulation engine; took {time.time() - start} secs.\n'
        min_logger.info(logger_string)
    except Exception as e:
        traceback.print_exc()
        min_logger.exception(f"Preparation of the simulation engine failed. Killing this worker's task.")
        clean_logger(min_logger)
        os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        
        return platform.node(), worker.id, full_start_time, time.time(), f'{pdb_file} failed on simulation preparation.'

    # run the minimization protocol and output minimized structure
    try:
        start = time.time()
        final_pdb, logger_string = run_minimization(simulation, path + model_descriptor)
        logger_string += f'Finished running the minimization calculation; took {time.time() - start} secs.\n'
        min_logger.info(logger_string)
    except Exception as e:
        traceback.print_exc()
        min_logger.exception(f"Simulation failed. Killing this worker's task.")
        clean_logger(min_logger)
        os.chmod(path + model_descriptor + '_min.log', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        
        return platform.node(), worker.id, full_start_time, time.time(), f'{pdb_file} failed on simulation.'
    
    del simulation

    return platform.node(), worker.id, full_start_time, time.time(), final_pdb
    






