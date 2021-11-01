#!/usr/bin/env python3
"""
    Creating the model minimization process used in AF to be used for post-modeling AF runs
    
"""
#######################################
### PREAMBLE
#######################################
import sys

from alphafold.common import protein
from alphafold.relax  import relax
from alphafold.relax  import utils
import jax.numpy as jnp

#######################################
### FUNCTIONS
#######################################


def pdb_to_string(pdb_file):
    """open and read pdb file, grabbing atom lines
    :param pdb_file: string, global or local path to a pdb file 
    :return: string of atomic information pulled from the pdb file
    """
    lines = []
    for line in open(pdb_file,"r"):
        if line[:6] == "HETATM" and line[17:20] == "MSE":
            line = "ATOM  "+line[6:17]+"MET"+line[20:]
        if line[:4] == "ATOM":
            lines.append(line)
    return "".join(lines)

#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # point to the file to be minimized
    example = sys.argv[1]
    chain_id = sys.argv[2]  # usually just 'A'
    # read in the pdb file to be minimized
    protein_obj = protein.from_pdb_string(pdb_to_string(example),chain_id=chain_id)
    # prepping the amber_relaxer method to be used.
    amber_relaxer = relax.AmberRelaxation(
              max_iterations=0,
              tolerance=2.39,
              stiffness=10.0,
              exclude_residues=[],
              max_outer_iterations=20)
    # minimize the protein
    relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=protein_obj)
    # write out the minimized structure to a pdb file
    with open("tmp.pdb", 'w') as f:
        f.write(relaxed_pdb_lines)




