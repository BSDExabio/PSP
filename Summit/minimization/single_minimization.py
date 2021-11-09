#!/usr/bin/env python3
"""
    Creating the model minimization process used in AF to be used for post-modeling AF runs
"""

#######################################
### PREAMBLE
#######################################
import sys
import time
import logging
import pdbfixer
from simtk import openmm
from simtk import unit
from simtk.openmm import app as openmm_app
from simtk.openmm.app.internal.pdbstructure import PdbStructure


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


def prep_protein(pdb_file,restraint_set,forcefield="amber99sb.xml",restraint_stiffness=10.0,exclude_residues=[],platform='CUDA',energy_units=unit.kilocalories_per_mole,length_units=unit.angstroms):
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
    
    logging.info(f"{delta_e} kcal mol^{-1} energy change ({einit},{efinal})\n")
    
    if not minimized:
        raise ValueError(f"Minimization failed after {fail_attempts} attempts.")
    

def run_pipeline(pdb_file,restraint_set,relax_exclude_residues):
    """
    """


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments
    pdb_file = sys.argv[1]
    restraint_set = sys.argv[2]
    relax_exclude_residues = sys.argv[3]

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


        

