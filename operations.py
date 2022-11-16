import numpy as np
"""

this script performs everyday operations, like conversion from 1D to 3D 
coordinates with periodic boundary conditions
@author: mrahm32
"""

def index_conversion_1D_3D(N, arrIndex):
    L = N**(1/3.0)
    L = np.round(L)
    i = arrIndex//(L**2)
    j = (arrIndex - i*(L**2))//L
    k = arrIndex - (i*(L**2) + j*L)
    coord = (i,j,k)
    return coord

def index_conversion_3D_1D(N,coord):
    
    L = N**(1/3.0)
    L = np.round(L)

    i, j, k = coord[0], coord[1], coord[2]
    arrIndex = i*(L**2) + j*(L) + k

    i_p, j_p, k_p = index_conversion_1D_3D(N, arrIndex)
    assert (i_p, j_p, k_p) == (i, j, k)
    return arrIndex

def int_2_str(matrix):
    str_matrix = []
    for row in matrix:
        new_row = np.array([int(row[0]), int(row[1]), int(row[2])])
        str_matrix.append(str(new_row))
    return str_matrix

def str_2_int(string_val):
    coord = string_val.split('[')[1].split(']')[0].split(' ')
    coord = np.array(coord)
    mask_space = (coord != '')
    coord = coord[mask_space]
    coord = [np.int(coord[0]), np.int(coord[1]), np.int(coord[2])]
    return coord
