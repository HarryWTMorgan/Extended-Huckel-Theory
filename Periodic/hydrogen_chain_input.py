"""
Extended Huckel Theory
Input file for a 1D hydrogen chain
"""

import numpy as np
import re

"""
First, create the unit cell in the same way as for the molecular calculations

Note - this way of naming and accessing atoms may be inefficient for generating
supercells. Consider reworking if you get time, but for now don't fix it until
it breaks!
"""

def create_H2_unit_cell():
    # Create dictionary of atoms and numerical indices
    
    atoms = {
    # add each atom-index pair to the dictionary
    "H1":0,
    "H2":1
    }
    
    # Create matrix of positions
    
    positions = np.array([[0.0,0.0,0.0], [0.0,0.0,0.74]])
    
    # Check for errors in positions/atoms input
    
    if len(atoms) < len(positions):
        print("Error: More atoms have been entered than coordinates")
    elif len(atoms) > len(positions):
        print("Error: More coordinates have been entered than atoms")
    
    # Put atomic orbital basis functions on each atom
    # array elements: basis function number, atom, orbital, n, l, m, zeta, ionization energy
    
    basis = np.array([["1","H1","1s",1,0,0,1.200,-0.5],\
                     ["2","H2","1s",1,0,0,1.200,-0.5]])

        
    return atoms,positions,basis



"""
We need to generate a supercell in order to calculate the overlaps between
the 'origin' unit cell and its nearest neighbour. We therefore need to extend
our list of atoms and coordinates.

1) List of atoms
    The atoms dictionary here is indexed numerically and by element; since
    element names can have different lengths we separate the two parts
    with a regular expression so we can add to the numbers while simply
    copying the element

2) Coordinates
    Much simpler than the list of atoms, especially in 1D.
    Loop through the unit cell positions and add a new position plus the
    lattice parameter to each one.
    
3) Basis
    To minimise the risk of errors creeping in at this stage, multiply the
    basis so that every atom has its own entry.
    If we had a different system of indexing atoms it might be possible to
    avoid expanding the basis (or even including each element only once).
"""

def generate_1D_supercell(atoms,positions,basis,a):
    """
    Generate new list of atoms
    Separate the label "H1" into a letter and a number, stored in a list.
    Do this for every atom label in the unit cell and store the resulting lists
    in split_atoms.
    Loop through the original atoms, duplicating the element and adding to the
    index.
    Currently the index of the new atom is the original plus the total
    number of atoms - this works for unit cells of one element
    (e.g. H1, H2 becomes H1, H2, H3, H4) but does not work so well for
    multiple elements (e.g. H1, He1 becomes H1, He1, H3, He3).
    Fixing this would require counting each element in the basis.
    This problem shouldn't cause any issues because each atom will still
    have a unique index
    """    
    split_atoms = []
    new_atoms = {}
    for atom in atoms:
        index = atom
        letter_and_number = re.split('(\d+)',index)
        del letter_and_number[-1]
        split_atoms.append(letter_and_number)
    for atom in range(len(atoms)):
        element = split_atoms[atom][0]
        number = int(split_atoms[atom][1])
        split_atoms.append([element,number+len(atoms)])
    for atom in range(len(split_atoms)):
        new_label = str(str(split_atoms[atom][0])+str(split_atoms[atom][1]))
        new_atoms[new_label] = atom
    print(new_atoms)

    """
    Generate new list of positions
    Loop through previous list, adding new coordinates with the lattice
    parameter added on
    For 1D case, assume periodicity along z
    """
    new_positions = positions
    for atom in range(len(positions)):
        new_positions = np.append(new_positions,[[positions[atom,0],positions[atom,1],positions[atom,2]+a]],axis=0)
    print(new_positions)
    
    """
    Generate new basis
    Double original basis set, then loop through the new basis modifying the
    basis function numbers and atom labels
    """
    new_basis = np.append(basis,basis,axis=0)
    atom_list = []
    for atom in new_atoms:
        atom_list.append(atom)
    for entry in range(len(new_basis)):
        new_basis[entry][0] = int(int(entry) + 1)
        new_basis[entry][1] = atom_list[entry]
    print(new_basis)
        
    




