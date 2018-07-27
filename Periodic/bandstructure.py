"""
Extended Huckel Theory
Calculate MOs and energy levels by diagonalising the Hamiltonian matrix
"""

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

from H2_chain_input import create_unit_cell, generate_1D_supercell
from periodic_overlap_matrix import k_overlap_matrix,construct_overlap_matrix
from periodic_hamiltonian_matrix import k_hamiltonian_matrix

UC_atoms, UC_positions, UC_basis = create_unit_cell()

a = 3.2

atoms, positions, basis = generate_1D_supercell(UC_atoms,UC_positions,UC_basis,a)


def band_structure(atoms,positions,basis,a,n):
    # Generate n uniformly distributed values in the range -pi/a to pi/a
    k_points = np.linspace((np.pi * -2/a),(np.pi * 2/a),n)
    band_energies = []
    MO_coeffs = []
    # Get overlap matrix for supercell 
    full_overlap_matrix = construct_overlap_matrix(atoms,positions,basis)
#    print(full_overlap_matrix)
    # Calculate energy levels for each k point
    for point in k_points:
        overlap_matrix = k_overlap_matrix(atoms, positions, basis, a, point)
        hamiltonian = k_hamiltonian_matrix(full_overlap_matrix, basis, a, point)
        rounded_MOs, energies = ordered_MOs(hamiltonian,overlap_matrix)
        eV_energies = energies_in_eV(energies)
        band_energies.append(eV_energies)
        MO_coeffs.append(rounded_MOs)
    return k_points, band_energies, MO_coeffs

def ordered_MOs(H,S):
    unsorted_energies, unsorted_MOs = scipy.linalg.eigh(H,S)
    indices = np.argsort(unsorted_energies)
    energies = unsorted_energies[indices]
    MOs = unsorted_MOs[:,indices]
    rounded_MOs = np.around(MOs, decimals=1)
    return rounded_MOs, energies
    

#MOs, energies = ordered_MOs(hamiltonian,overlap_matrix)

# Convert energies to eV and round

def energies_in_eV(energy_list):
    energy_list = np.real(energy_list)
    for E in energy_list:
        E *= 27.2114
    energy_list = np.around(energy_list, decimals=5)
    return energy_list

#eV_energies = energies_in_eV(energies)


k_points,band,MOs = band_structure(atoms,positions,basis,a,200)

plt.plot(k_points,band)
"""
Use Giacomo Marchioro's energy level diagram script to plot MO diagram


from energydiagram import ED
diagram = ED()
for level in eV_energies:
    diagram.add_level(level,'' ,'last')
diagram.plot()
"""