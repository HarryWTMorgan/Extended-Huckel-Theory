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
    S_ij = complex(overlaps[i,j])
    H_ii = complex(basis[i,-1])
    H_jj = complex(basis[j,-1])
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
    hamiltonian = np.zeros((len(overlaps),len(overlaps)),dtype=complex)
    for i in range(len(overlaps)):
        for j in range(i+1): # Compute lower triangle
            # Diagonal elements are on-site ionization energies
            if i == j :
                hamiltonian[i,j] = float(basis[i,-1])
            # Off-diagonal elements given by Wolfsberg-Helmholtz formula
            else:
                hamiltonian[i,j] = wolfsberg_helmholtz(overlaps, basis, i, j)
    # Reflect in main diagonal to get symmetric matrix
            hamiltonian[j,i] = hamiltonian[i,j]
    return hamiltonian


def slice_hamiltonian_matrix(overlaps,basis):
    full_hamiltonian_matrix = compute_hamiltonian_matrix(overlaps,basis)
    H0 = full_hamiltonian_matrix[:int(len(basis)/2),:int(len(basis)/2)]
    H1 = full_hamiltonian_matrix[:int(len(basis)/2),int(len(basis)/2):]
    H2 = full_hamiltonian_matrix[int(len(basis)/2):,:int(len(basis)/2)]
    return H0,H1,H2

        
def k_hamiltonian_matrix(overlaps,basis,a,k):
    # Get the useful parts of the hamiltonian matrix
    H0,H1,H2 = slice_hamiltonian_matrix(overlaps,basis)
    # Turn H0 into a complex array (even though it's real)
    new_H = np.zeros((len(H0),len(H0)),dtype=complex)
    new_H += H0
    
    pos_phase_factor = np.exp(1j*k*a)
    neg_phase_factor = np.exp(-1j*k*a)
#    print("H0",new_H,'\n',"H1",H1,'\n',"H2",H2)
    new_H = new_H + H1 * pos_phase_factor + H2 * neg_phase_factor
#    print("new H",new_H,'\n')
    return(new_H)
     