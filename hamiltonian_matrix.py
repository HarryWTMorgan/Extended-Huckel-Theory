"""
Extended Huckel Theory
Compute Hamilton matrix from overlap matrix
"""

import numpy as np

"""
Wolfsberg-Helmholtz formula for off-diagonal Hamiltonian matrix elements
Function takes in overlap matrix S, ionisation energies and atomic orbital indices i,j
"""

def wolfsberg_helmholtz(overlaps,basis,i,j):
    K = 1.75
    S_ij = float(overlaps[i,j])
    H_ii = float(basis[i,-1])
    H_jj = float(basis[j,-1])
    H_ij = K * S_ij * 0.5 * (H_ii + H_jj)
    return H_ij

"""
Use overlap matrix for methane from Lowe as test case
"""

CH4_overlaps = np.array([[1.0, 0.0, 0.0, 0.0, 0.5133, 0.5133, 0.5133, 0.5133],\
                        [0.0, 1.0, 0.0, 0.0, 0.4855, -0.1618, -0.1618, -0.1618],\
                        [0.0, 0.0, 1.0, 0.0, 0.0, 0.4577, -0.2289, -0.2289],\
                        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.3964, -0.3964],\
                        [0.5133, 0.4855, 0.0, 0.0, 1.0, 0.1805, 0.1805, 0.1805],\
                        [0.5133, -0.1618, 0.4577, 0.0, 0.1805, 1.0, 0.1805, 0.1805],\
                        [0.5133, -0.1618, -0.2289, 0.3964, 0.1805, 0.1805, 1.0, 0.1805],\
                        [0.5133, -0.1618, -0.2289, -0.3964, 0.1805, 0.1805, 0.1805, 1.0]])


"""
Generate Hamiltonian matrix from overlap matrix
Ionization energies in basis set
"""

def compute_hamiltonian_matrix(overlaps, basis):
    # Create matrix of zeros of the appropriate size
    hamiltonian = np.zeros((len(overlaps),len(overlaps)))
    for i in range(len(overlaps)):
        for j in range(i+1): # Compute lower triangle
            # Diagonal elements are on-site ionization energies
            if i == j :
                hamiltonian[i,j] = basis[i,-1]
            # Off-diagonal elements given by Wolfsberg-Helmholtz formula
            else:
                hamiltonian[i,j] = wolfsberg_helmholtz(overlaps, basis, i, j)
    # Reflect in main diagonal to get symmetric matrix
            hamiltonian[j,i] = hamiltonian[i,j]
    return hamiltonian
 