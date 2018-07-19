"""
Extended Huckel Theory
Compute overlap matrix
"""

#from input_molecule import create_methane

import numpy as np
from math import factorial as bang

#atoms,positions,basis = create_ammonia()

"""
First we have some short functions to calculate particular quantities
in the radial overlap integral below, written seperately for simplicity.
Functions to calculate:
    Scalar distance between two atoms
    Numerical factor before all the sums in the radial overlap expression
    C_lmj coefficients describing radial overlap between different orbital types
    A and B terms (from Mulliken et al)
"""

# Calculate interatomic distance between atoms a and b in atomic units

def bond_length(positions, a, b):
    # Conversion for Angstrom to bohr (0.529 bohr per Angstrom)
    A0=0.52917721092
    bond_vector = positions[b] - positions[a]
    length = np.linalg.norm(bond_vector)
    length = length / A0
    return length

# Calculate pre-factor in the radial overlap equation
# (the first term, outside the sum)
# for atoms a, b

def calc_prefactor(atoms,positions,basis,a,b):
    # Extract orbital properties from the basis set
    zeta_a = float(basis[a,-2])
    zeta_b = float(basis[b,-2])
    n_a = float(basis[a,3])
    n_b = float(basis[b,3])
    l_a = float(basis[a,4])
    l_b = float(basis[b,4])
    
    # Use atoms dictionary to match orbitals to atoms,
    # then get their coordinates from the list of positions
    # to calculate the bond length
    # label = string, index = integer
    
    atom_label_a = basis[a,1]
    atom_label_b = basis[b,1]
    atom_index_a = atoms[atom_label_a]
    atom_index_b = atoms[atom_label_b]
    R_ab = bond_length(positions,atom_index_a,atom_index_b)
    # bang = factorial
    prefactor = 0.5 * (zeta_a ** (n_a + 1/2) * (zeta_b ** (n_b + 1/2))) * \
                ((((2*l_a + 1) * (2*l_b + 1)) / (bang(2*n_a) * bang(2*n_b))) ** 0.5) \
                * R_ab ** (n_a + n_b + 1)
    return prefactor

# Write and access a dictionary of C_lmj coefficients
# Add more values to the dictionary if you want to use d orbitals    

def get_C_lmj(l,m,j):
    # Create dictionary
    C_lmj_dict = {
                (0,0,0) : 1,
                (1,0,0) : 1,
                (1,1,0) : 0.5**0.5
                }
    # Get the relevant value from the dictionary
    C_lmj = C_lmj_dict[(l,m,j)]
    return C_lmj

# Formulae for A and B quantities found in radial overlap expression
# Taken from Mulliken, Rieke, Orloff & Orloff, J. Chem. Phys 17, 1248 (1949)

def A_func(k,p):
    # Initialise A at zero
    A = 0
    # Sum over mu - upper limit is k+2 so the expression is evaluated for mu=k+1
    for mu in range(1, k+2):
        A += bang(k) / ((p ** mu) * bang(k - mu + 1))
    # Now multiply by the exponential term
    A *= np.exp(-p)
    return A

# Argument is written as 'pt' here so that the formula looks like that given
# by Mulliken - when we come to put this function into the radial overlap
# equation, pt will be replaced by p * t
    
def B_func(k,pt):
    # Initialise B at zero
    B = 0
    # B(0) can be evaluated much more straightforwardly (see Mulliken or Stewart)
    # so include a check here to avoid unnecessary calculations
    if pt == 0:
        if k % 2 == 1: # k odd
            B = 0
        elif k % 2 == 0: # k even
            B = 2 / (k+1)
    else:
        # Calculate each of the sums seperately for ease of understanding, then
        # multiply them by their exponential factors, then subtract them from B
    
        # Label the two sums B_1 and B_2 - calculate them both in the same loop
        # since the sums run over the same indices
        
        B_1 = 0
        B_2 = 0
        for mu in range(1, k+2):
            B_1 += bang(k) / ((pt ** mu) * bang(k - mu + 1))
            B_2 += ( bang(k) * (-1)**(k-mu) ) / ((pt ** mu) * bang(k - mu + 1))
        B -= ( B_1 * np.exp(-pt) + B_2 * np.exp(pt) )
    return B

"""
Now we write the expression for the radial overlap between two orbitals.
This consists of sums over 8 indices - 4 for each atom (a and b)

This formula comes from the article "Integrals: Overlap" by James J. P. Stewart
in the Encyclopedia of Computational Chemistry, vol. 2

m (magnetic quantum number) is given to the function directly rather than
extracted from the basis set - this is because for the purposes of this
calcuation we align the orbitals along their own principal axes so it
is not practical to use the m values in the basis set.
Working this way enables us to calculate overlaps between s and px/py orbitals,
and sigma/pi overlaps for a pair of p orbitals.
"""
    
def calc_radial_overlap(atoms, positions,basis,a,b,m):
    # Extract orbital properties from the basis set
    zeta_a = float(basis[a,-2])
    zeta_b = float(basis[b,-2])
    na = int(basis[a,3])
    nb = int(basis[b,3])
    la = int(basis[a,4])
    lb = int(basis[b,4])
    
    # Use atoms dictionary to match orbitals to atoms,
    # then get their coordinates from the list of positions
    # to calculate the bond length
    # label = string, index = integer
    
    atom_label_a = basis[a,1]
    atom_label_b = basis[b,1]
    atom_index_a = atoms[atom_label_a]
    atom_index_b = atoms[atom_label_b]
    R_ab = bond_length(positions,atom_index_a,atom_index_b)
    
    # Define p and t
    
    p = (R_ab * (zeta_a + zeta_b)) / 2
    t = (zeta_a - zeta_b) / (zeta_a + zeta_b)
    
    # Initialise overlap at 0
    radial_overlap = 0
    # Write each sum as a for loop
    for ja in range(int((la-m)/2)+1):
        for jb in range(int((lb-m)/2)+1):
            for ka in range(m+1):
                for kb in range(m+1):
                    for Pa in range(na-la-2*ja+1):
                        for Pb in range(nb-lb-2*jb+1):
                            for qa in range(la-m-2*ja+1):
                                for qb in range(lb-m-2*jb+1):
                                    # Write each term on a different line and
                                    # multiply them all together
                                    radial_overlap += \
                                    get_C_lmj(la,m,ja)\
                                    * \
                                    get_C_lmj(lb,m,jb)\
                                    * \
                                    ((bang(lb-m-2*jb)) / (bang(lb-m-2*jb-qb)*bang(qb)))\
                                    * \
                                    ((bang(la-m-2*ja)) / (bang(la-m-2*ja-qa)*bang(qa)))\
                                    * \
                                    ((bang(na-la+2*ja)) / (bang(na-la+2*ja-Pa)*bang(Pa)))\
                                    * \
                                    ((bang(nb-lb+2*jb)) / (bang(nb-lb+2*jb-Pb)*bang(Pb)))\
                                    * \
                                    ((bang(m))**2) / (bang(m-ka)*bang(ka)*bang(m-kb)*bang(kb))\
                                    * \
                                    (-1) ** (ka+kb+m+Pb+qb) \
                                    * \
                                    B_func(2*ka+Pa+Pb+qa+qb,p*t) \
                                    * \
                                    A_func(2*kb+na-la+2*ja+nb-lb+2*jb-Pa-Pb+qa+qb,p)
    radial_overlap *= calc_prefactor(atoms, positions,basis,a,b)
    return radial_overlap                                    


"""
The next step is to rotate the orbitals, which have so far been perfectly aligned,
into their molecular orientation. This involves calculating the "direction cosines"
(projections of the bond vector along the coordinate axes) and looking up the
appropriate formula from a list according to the orbital symmetries
"""

# Calculate direction cosines
# (Cartesian components of normalised vector)
def calc_direction_cosines(atoms, positions, basis, a, b):
    # Use atoms dictionary to match orbitals to atoms,
    # then get their coordinates from the list of positions
    # to calculate the bond length
    # label = string, index = integer
    
    atom_label_a = basis[a,1]
    atom_label_b = basis[b,1]
    atom_index_a = atoms[atom_label_a]
    atom_index_b = atoms[atom_label_b]
    
    bond_vector = positions[atom_index_b] - positions[atom_index_a]
    bond_length = np.linalg.norm(bond_vector)
    unit_bond_vec = bond_vector / bond_length
    l,m,n = unit_bond_vec[0],unit_bond_vec[1],unit_bond_vec[2]
    return l,m,n

"""
This function looks up the orbital quantum numbers in the basis set and uses
them to select the appropriate formula

Taken from https://en.wikipedia.org/wiki/Tight_binding#Table_of_interatomic_matrix_elements
"""
    
def calc_overlap(atoms, positions, basis, a, b):
    # Extract orbital properties from the basis set
    la = int(basis[a,4])
    lb = int(basis[b,4])
    ma = int(basis[a,5])
    mb = int(basis[b,5])
    
    # Calculate Cartesian components of R_ab
    l,m,n = calc_direction_cosines(atoms, positions, basis, a, b)

    
    # ss overlap
    if la == 0 and lb == 0:
        overlap = calc_radial_overlap(atoms, positions, basis, a, b, 0)
        # no angular dependence
    
    # sp and ps overlap - ps is +ve and sp is -ve since pz is
    # defined as +ve in the +z direction
    elif (la == 0 and lb == 1) or (la == 1 and lb == 0):
        # Calculate radial part
        radial_overlap = calc_radial_overlap(atoms, positions, basis, a, b, 0)
        
        # Find which atom the p orbital is on
        # ps
        if la == 1:
            mp = ma

        # sp
        elif lb == 1:
            mp = mb
            # change sign of overlap
            radial_overlap = -radial_overlap
        
        # project on to the appropriate axis                
        # s px
        if mp == 1:
            overlap = l * radial_overlap
        # s py
        elif mp == -1:
            overlap = m * radial_overlap
        # s pz
        elif mp == 0:
            overlap = n * radial_overlap
    
    # pp overlap
    elif la == 1 and lb == 1:
        # Calculate radial parts for sigma and pi symmetry
        # NB sigma overlap is -ve
        sigma_overlap = -calc_radial_overlap(atoms, positions, basis, a, b, 0)
        pi_overlap = calc_radial_overlap(atoms, positions, basis, a, b, 1)
        
        # px on a
        if ma == 1:
            # px px
            if mb == 1:
                overlap = l**2 * sigma_overlap + (1 - l**2) * pi_overlap
            # px py
            if mb == -1:
                overlap = l*m * sigma_overlap + (1 - l*m) * pi_overlap
            # px pz
            if mb == 0:
                overlap = l*n * sigma_overlap + (1 - l*n) * pi_overlap
        
        # py on a
        elif ma == -1:
            # py px
            if mb == 1:
                overlap = l*m * sigma_overlap + (1 - l*m) * pi_overlap           
            # py py
            if mb == -1:
                overlap = m**2 * sigma_overlap + (1 - m**2) * pi_overlap
            # py pz
            if mb == 0:
                overlap = m*n * sigma_overlap + (1 - m*n) * pi_overlap
        
        # pz on a
        elif ma == 0:
            # pz px
            if mb == 1:
                overlap = l*n * sigma_overlap + (1 - l*n) * pi_overlap           
            # pz py
            if mb == -1:
                overlap = m*n * sigma_overlap + (1 - m*n) * pi_overlap
            # pz pz
            if mb == 0:
                overlap = n**2 * sigma_overlap + (1 - n**2) * pi_overlap
                        
    return overlap
      
"""
Finally, generate the overlap matrix - each element is an overlap between two
orbitals in the basis. Here we create a square zero matrix of the right size,
set the diagonal elements to 1, set any overlaps between orbitals on the same
atom to zero (enforcing orthogonality) and then evaluate all other elements.
"""


def construct_overlap_matrix(atoms, positions,basis):
    # Create square matrix of zeros of the same length as the basis set
    overlap_matrix = np.zeros((len(basis),len(basis)))
    # Loop through all rows
    for i in range(len(overlap_matrix)):
        # Compute lower triangle only
        for j in range(i+1):
            # Diagonal elements are 1 (orbitals fully overlap with themselves!)
            if i == j:
                overlap_matrix[i,j] = 1.0
            # Enforce orthogonality of atomic orbitals on the same atom
            elif basis[i,1] == basis[j,1]:                
                overlap_matrix[i,j] = 0.0
            # Orbitals on different atoms overlap according to formula TBA!
            else:
                overlap_matrix[i,j] = calc_overlap(atoms, positions, basis, i, j)
            # Symmetrize the matrix
            overlap_matrix[j,i] = overlap_matrix[i,j]
    return overlap_matrix

#overlaps = construct_overlap_matrix(positions,basis)

