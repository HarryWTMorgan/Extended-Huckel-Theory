"""
Extended Huckel Theory
Combine input, overlap, hamiltonian and MO scripts to give MOs and energies
"""

from input_molecule import create_methane
from overlap_matrix import *
from hamiltonian_matrix import wolfsberg_helmholtz, compute_hamiltonian_matrix

atoms, positions, basis = create_methane()

overlap_matrix = construct_overlap_matrix(positions, basis)

hamiltonian_matrix = compute_hamiltonian_matrix