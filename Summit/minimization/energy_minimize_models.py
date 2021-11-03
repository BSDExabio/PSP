#!/usr/bin/env python3
"""
    Creating the model minimization process used in AF to be used for post-modeling AF runs
    
"""

#######################################
### PREAMBLE
#######################################
import sys
import time
import numpy as np
from simtk import openmm
from simtk import unit
from simtk.openmm import app as openmm_app
from simtk.openmm.app.internal.pdbstructure import PdbStructure

ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms
RELAX_MAX_ITERATIONS = 0        # fed to openmm minimizeEnergy
RELAX_ENERGY_TOLERANCE = 2.39   # fed to openmm minimizeEnergy
RELAX_STIFFNESS = 10.0          # fed to openmm minimizeEnergy
RELAX_EXCLUDE_RESIDUES = []     # fed to openmm minimizeEnergy; left empty by default
MAX_FAIL_ATTEMPTS = 100         # used to limit number of attempts to successfully run minimization
RELAX_MAX_OUTER_ITERATIONS = 20 # used to cutoff outer while 

tolerance = tolerance * ENERGY  # 2.39 kcal mol^{-1}
stiffness = 10  * ENERGY / (LENGTH**2)  # 10 kcal mol^{-1} \AA^{-2}
restraint_set = "non_hydrogen"      # or "c_alpha"
exclude_residues = []               # list of residue indices

#######################################
### COMMAND LINE ARGUMENTS
#######################################
pdb_file = sys.argv[1]

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
    system: openmm.System,
    reference_pdb: openmm_app.PDBFile,
    stiffness: unit.Unit,
    rset: str,
    exclude_residues: Sequence[int]):
  """Adds a harmonic potential that restrains the end-to-end distance."""
  assert rset in ["non_hydrogen", "c_alpha"]

  force = openmm.CustomExternalForce(
      "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
  force.addGlobalParameter("k", stiffness)
  for p in ["x0", "y0", "z0"]:
    force.addPerParticleParameter(p)

  for i, atom in enumerate(reference_pdb.topology.atoms()):
    if atom.residue.index in exclude_residues:
      continue
    if will_restrain(atom, rset):
      force.addParticle(i, reference_pdb.positions[i])
  logging.info("Restraining %d / %d particles.",
               force.getNumParticles(), system.getNumParticles())
  system.addForce(force)

#######################################
### MAIN
#######################################

# load pdb file into an openmm Topology and coordinates object
pdb = openmm_app.PDBFile(pdb_file)
# set the FF and constraints objects
force_field = openmm_app.ForceField("amber99sb.xml")
constraints = openmm_app.HBonds
# create the simulation ready system
system = force_field.createSystem(pdb.topology, constraints=constraints)
# prep and add restraints to the simulation system
if stiffness > 0 * ENERGY / (LENGTH**2):
    # code up the _add_restraints function
    _add_restraints(system, pdb, stiffness, restraint_set, exclude_residues)

# required to set this for prepping the simulation object
integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
# determine what hardware will be used to perform the calculations
platform = openmm.Platform.getPlatformByName("CUDA")
#platform = openmm.Platform.getPlatformByName("CPU")
# prep the simulation object
simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)

outer_start = time.time()
# loop over max_outer_iterations
while iteration < RELAX_MAX_OUTER_ITERATIONS:   # ignoring violations, which may decrease efficicency of this pipeline relative to AF's

    # set the atom positions for the simulation's system's topology
    simulation.context.setPositions(pdb.positions)
        # need to check what shape/format that the pdb.positions has; maybe I can just read in the ret["min_pdb"] 

    inner_start = time.time()
    attempts = 0
    minimized = False
    while not minimized and attempts < MAX_FAIL_ATTEMPTS:
        attempts += 1
        try:
            # return energies and positions
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            einit = state.getPotentialEnergy().value_in_unit(ENERGY)
            posinit = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
            # running minimization
            simulation.minimizeEnergy(maxIterations=max_iterations,tolerance=tolerance)
            # return energies and positions
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            efinal = state.getPotentialEnergy().value_in_unit(ENERGY)
            pos = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
            min_pdb = _get_pdb_string(simulation.topology, state.getPositions())
        except Exception as e:
            print(e)
    inner_time_end = time.time() - inner_start
    delta_e = 
    print( f"{inner_time_end} secs spent on minimization {iteration}\n {delta_e} kcal mol^{-1} \AA^{-2} energy change\n")
    if not minimized:
        raise ValueError(f"Minimization failed after {MAX_FAIL_ATTEMPTS} attempts.")
    iteration += 1

# calculate rmsd between pos and posinit
# write out the min_pdb


