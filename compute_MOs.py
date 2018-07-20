"""
Extended Huckel Theory
Calculate MOs and energy levels by diagonalising the Hamiltonian matrix
"""

import numpy as np
import scipy.linalg
#import matplotlib.pyplot as plt

from dinitrogen_input import create_molecule
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
    rounded_MOs = np.around(MOs, decimals=1)
    return rounded_MOs, energies
    

MOs, energies = ordered_MOs(hamiltonian,overlap_matrix)

# Convert energies to eV and round

def energies_in_eV(energy_list):
    energy_list = np.real(energy_list)
    for E in energy_list:
        E *= 27.2114
    energy_list = np.around(energy_list, decimals=3)
    return energy_list

eV_energies = energies_in_eV(energies)

"""
Use Giacomo Marchioro's energy level script to plot MO diagram
"""

from energydiagram import ED
diagram = ED()
for level in eV_energies:
    diagram.add_level(level,'' ,'last')
diagram.plot()