"""
Extended Huckel Theory
Calculate MOs and energy levels by diagonalising the Hamiltonian matrix
"""

import numpy as np
import scipy.linalg
#import matplotlib.pyplot as plt

from ammonia_input import create_ammonia
from overlap_matrix import construct_overlap_matrix
from hamiltonian_matrix import compute_hamiltonian_matrix

atoms, positions, basis = create_ammonia()

overlap_matrix = construct_overlap_matrix(atoms, positions, basis)

hamiltonian = compute_hamiltonian_matrix(overlap_matrix, basis)

def ordered_MOs(H,S):
    unsorted_energies, unsorted_MOs = scipy.linalg.eig(H,S)
    indices = np.argsort(unsorted_energies)
    energies = unsorted_energies[indices]
    MOs = unsorted_MOs[:,indices]
    return MOs, energies
    

MOs, energies = ordered_MOs(hamiltonian,overlap_matrix)

