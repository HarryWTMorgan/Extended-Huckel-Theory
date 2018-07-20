"""
Extended Huckel Theory
Nitrogen molecule (N2) input
"""

import numpy as np


def create_molecule():
    # Create dictionary of atoms and numerical indices
    
    atoms = {
    # add each atom-index pair to the dictionary
    "N1":0,
    "N2":1
    }
    
    # Create matrix of positions
    
    positions = np.array([[0.0,0.0,0.0], [0.0,0.0,1.1]])
    
    # Check for errors in positions/atoms input
    
    if len(atoms) < len(positions):
        print("Error: More atoms have been entered than coordinates")
    elif len(atoms) > len(positions):
        print("Error: More coordinates have been entered than atoms")
    
    # Put atomic orbital basis functions on each atom
    # array elements: basis function number, atom, orbital, n, l, m, zeta, ionization energy
    
    basis = np.array([["1","N1","2s",2,0,0,2.140,-0.9554],\
                     ["2","N1","2pz",2,1,0,1.950,-0.4925],\
                     ["3","N1","2px",2,1,1,1.950,-0.4925],\
                     ["4","N1","2py",2,1,-1,1.950,-0.4925],\
                     ["5","N2","2s",2,0,0,2.140,-0.9554],\
                     ["6","N2","2pz",2,1,0,1.950,-0.4925],\
                     ["7","N2","2px",2,1,1,1.950,-0.4925],\
                     ["8","N2","2py",2,1,-1,1.950,-0.4925]])

        
    return atoms,positions,basis

