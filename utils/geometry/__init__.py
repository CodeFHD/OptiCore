"""
Copyright 2019-2024, Johannes Hinrichs

This file is part of OptiCore.

OptiCore is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OptiCore is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OptiCore. If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

"""
To apply rotation matrix R to a vector V of shape (3,):
np.matmul(R, V)

to apply to a multiple vector U of shape (N, 3)
np.matmul(R, U.T).T
"""

def get_rotmat_x(phi):
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    R = [[1,  0,     0    ],
         [0,  cphi,  -sphi],
         [0,  sphi,  cphi ]]
    return np.array(R)

def get_rotmat_y(phi):
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    R = [[cphi,   0,  sphi],
         [0,      1,  0   ],
         [-sphi,  0,  cphi]]
    return np.array(R)

def get_rotmat_z(phi):
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    R = [[cphi,  -sphi,  0],
         [sphi,  cphi,   0],
         [0,     0,      1]]
    return np.array(R)

def get_rotmat_axis(phi, a=[1,0,0]):
    a1, a2, a3 = a
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    mcphi = 1 - cphi
    R11 = a1*a1*mcphi + cphi
    R12 = a1*a2*mcphi - a3*sphi
    R13 = a1*a3*mcphi + a2*sphi
    R21 = a1*a2*mcphi + a3*sphi
    R22 = a2*a2*mcphi + cphi
    R23 = a2*a3*mcphi - a1*sphi
    R31 = a1*a3*mcphi - a2*sphi
    R32 = a2*a3*mcphi + a1*sphi
    R33 = a3*a3*mcphi + cphi
    R = [[R11, R12, R13],
         [R21, R22, R23],
         [R31, R32, R33]]
    return np.array(R)

def rotate_vector_x(V, surf_rotation):
    input_is_array = isinstance(V, np.ndarray)
    R_z = get_rotmat_x(surf_rotation)
    if input_is_array:
        V = np.matmul(R_z, V.T).T
    else:
        V = np.array(V)
        V = np.matmul(R_z, V)
        V = [V[0], V[1], V[2]]
    return V

def rotate_vector_y(V, surf_rotation):
    input_is_array = isinstance(V, np.ndarray)
    R_z = get_rotmat_y(surf_rotation)
    if input_is_array:
        V = np.matmul(R_z, V.T).T
    else:
        V = np.array(V)
        V = np.matmul(R_z, V)
        V = [V[0], V[1], V[2]]
    return V

def rotate_vector_z(V, surf_rotation):
    input_is_array = isinstance(V, np.ndarray)
    R_z = get_rotmat_z(surf_rotation)
    if input_is_array:
        V = np.matmul(R_z, V.T).T
    else:
        V = np.array(V)
        V = np.matmul(R_z, V)
        V = [V[0], V[1], V[2]]
    return V


""" TEST """

if __name__ == '__main__':
    print()
    testvector = np.arange(3)
    print('Testvector:')
    print(testvector)
    print()
    R = get_rotmat_z(np.pi/4)
    print('Rotmat:')
    print(R)
    print()
    result = np.matmul(R, testvector)
    print('Result:')
    print(result)
    print()
    print('-'*60)
    print()
    testvector = np.arange(5*3).reshape((5,3))
    print('Testvector (long):')
    print(testvector)
    print()
    R = get_rotmat_z(np.pi/4)
    print('Rotmat:')
    print(R)
    print()
    result = np.matmul(R, testvector.T).T
    print('Result:')
    print(result)
    print()
