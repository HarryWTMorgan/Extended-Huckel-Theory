"""
Extended Huckel theory
Input file for NH3
"""

import numpy as np


def create_molecule():
    # Create dictionary of atoms and numerical indices
    
    atoms = {
    # add each atom-index pair to the dictionary
    "N":0,
    "H1":1,
    "H2":2,
    "H3":3,
    }
    
    # Create matrix of positions
    
    positions = np.array([[0.257,-0.363,0.0], [0.257,0.727,0.0], [0.771,-0.727,0.890], \
                         [0.771, -0.727, -0.890]])
    
    # Check for errors in positions/atoms input
    
    if len(atoms) < len(positions):
        print("Error: More atoms have been entered than coordinates")
    elif len(atoms) > len(positions):
        print("Error: More coordinates have been entered than atoms")
    
    # Put atomic orbital basis functions on each atom
    # array elements: basis function number, atom, orbital, n, l, m, zeta, ionization energy
    
    basis = np.array([["1","N","2s",2,0,0,2.140,-0.9554],\
                     ["2","N","2pz",2,1,0,1.950,-0.4925],\
                     ["3","N","2px",2,1,1,1.950,-0.4924],\
                     ["4","N","2py",2,1,-1,1.950,-0.4924],\
                     ["5","H1","1s",1,0,0,1.200,-0.5],\
                     ["6","H2","1s",1,0,0,1.200,-0.5],\
                     ["7","H3","1s",1,0,0,1.200,-0.5]])

        
    return atoms,positions,basis

 