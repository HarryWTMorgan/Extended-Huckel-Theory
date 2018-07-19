"""
Extended Huckel Theory
Calculate MOs and energy levels by diagonalising the Hamiltonian matrix
"""

import numpy as np
import scipy.linalg
#import matplotlib.pyplot as plt

from ethene_input import create_molecule
from overlap_matrix import construct_overlap_matrix
from hamiltonian_matrix import compute_hamiltonian_matrix

atoms, positions, basis = create_molecule()

overlap_matrix = construct_overlap_matrix(atoms, positions, basis)

hamiltonian = compute_hamiltonian_matrix(overlap_matrix, basis)

def ordered_MOs(H,S):
    unsorted_energies, unsorted_MOs = scipy.linalg.eig(H,S)
    indices = np.argsort(unsorted_energies)
    energies = unsorted_energies[indices]
    MOs = unsorted_MOs[:,indices]
    rounded_MOs = np.around(MOs, decimals=2)
    return rounded_MOs, energies
    

MOs, energies = ordered_MOs(hamiltonian,overlap_matrix)

