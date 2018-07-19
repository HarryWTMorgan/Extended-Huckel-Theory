"""
Extended Huckel Theory
Calculate MOs and energy levels by diagonalising the Hamiltonian matrix
"""

import numpy as np
import scipy.linalg
#import matplotlib.pyplot as plt

from input_molecule import create_methane
from overlap_matrix import construct_overlap_matrix
from hamiltonian_matrix import compute_hamiltonian_matrix

atoms, positions, basis = create_methane()

overlap_matrix = construct_overlap_matrix(atoms, positions, basis)

hamiltonian = compute_hamiltonian_matrix(overlap_matrix, basis)

energies, MOs = scipy.linalg.eig(hamiltonian,overlap_matrix)

ordered_energies = np.sort(energies)

print(ordered_energies)

