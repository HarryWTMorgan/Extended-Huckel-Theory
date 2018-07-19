"""
Extended Huckel Theory
Ethene input file
"""

import numpy as np


def create_molecule():
    # Create dictionary of atoms and numerical indices
    
    atoms = {
    # add each atom-index pair to the dictionary
    "C1":0,
    "C2":1,
    "H1":2,
    "H2":3,
    "H3":4,
    "H4":5
    }
    
    # Create matrix of positions
    
    positions = np.array([[0.0,0.0,0.667], [0.0,0.0,-0.667], [0.0,0.923,1.238], \
                         [0.0,-0.923,1.238],[0.0,0.923,-1.238],[0.0,-0.923,-1.238]])
    
    # Check for errors in positions/atoms input
    
    if len(atoms) < len(positions):
        print("Error: More atoms have been entered than coordinates")
    elif len(atoms) > len(positions):
        print("Error: More coordinates have been entered than atoms")
    
    # Put atomic orbital basis functions on each atom
    # array elements: basis function number, atom, orbital, n, l, m, zeta, ionization energy
    
    basis = np.array([["1","C1","2s",2,0,0,1.625,-0.7144],\
                     ["2","C1","2pz",2,1,0,1.625,-0.3921],\
                     ["3","C1","2px",2,1,1,1.625,-0.3921],\
                     ["4","C1","2py",2,1,-1,1.625,-0.3921],\
                     ["5","C2","2s",2,0,0,1.625,-0.7144],\
                     ["6","C2","2pz",2,1,0,1.625,-0.3921],\
                     ["7","C2","2px",2,1,1,1.625,-0.3921],\
                     ["8","C2","2py",2,1,-1,1.625,-0.3921],\
                     ["9","H1","1s",1,0,0,1.200,-0.5],\
                     ["10","H2","1s",1,0,0,1.200,-0.5],\
                     ["11","H3","1s",1,0,0,1.200,-0.5],\
                     ["12","H4","1s",1,0,0,1.200,-0.5]])

        
    return atoms,positions,basis
