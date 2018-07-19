"""
Extended Huckel Theory
Input molecular coordinates and ionization energies
"""

import numpy as np


def create_methane():
    # Create dictionary of atoms and numerical indices
    
    atoms = {
    # add each atom-index pair to the dictionary
    "C":0,
    "H1":1,
    "H2":2,
    "H3":3,
    "H4":4
    }
    
    # Create matrix of positions
    
    positions = np.array([[0.0,0.0,0.0], [0.0,0.0,1.1], [1.03709,0.0,-0.366667], \
                         [-0.518545,0.898146,-0.366667],[-0.518545,-0.898146,-0.366667]])
    
    # Check for errors in positions/atoms input
    
    if len(atoms) < len(positions):
        print("Error: More atoms have been entered than coordinates")
    elif len(atoms) > len(positions):
        print("Error: More coordinates have been entered than atoms")
    
    # Put atomic orbital basis functions on each atom
    # array elements: basis function number, atom, orbital, n, l, m, zeta, ionization energy
    
    basis = np.array([["1","C","2s",2,0,0,1.625,-0.7144],\
                     ["2","C","2pz",2,1,0,1.625,-0.3921],\
                     ["3","C","2px",2,1,1,1.625,-0.3921],\
                     ["4","C","2py",2,1,-1,1.625,-0.3921],\
                     ["5","H1","1s",1,0,0,1.200,-0.5],\
                     ["6","H2","1s",1,0,0,1.200,-0.5],\
                     ["7","H3","1s",1,0,0,1.200,-0.5],\
                     ["8","H4","1s",1,0,0,1.200,-0.5]])

        
    return atoms,positions,basis

                     