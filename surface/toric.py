"""
Copyright 2019-2025, Johannes Hinrichs

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

from .surfaceutils import get_xy_rotatedsurface
from ..utils.geometry import rotate_vector_x, rotate_vector_z

def get_zN_toric(x, y, R, r, surf_rotation=0, coord_input='original', returnformat=''):
    if coord_input == 'original' and not surf_rotation == 0:
        # rotate the coordinates to align with the cylidner axis
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    # Here, the functions for z and N are combined for numerical efficiency
    R_neg = R < 0
    r_neg = r < 0
    sign_R = 1 - 2*R_neg
    sign_r = 1 - 2*r_neg
    sign_prod = sign_R*sign_r
    # sign_both = 1 - 2*(R_neg and r_neg)
    r = abs(r)
    R = abs(R)
    R_T = R - r*sign_prod

    # get z
    sqrt_1 = sign_prod*np.sqrt(r**2 - y**2)
    sqrt_2 = np.sqrt((R_T + sqrt_1)**2 - x**2)
    z = (R - sqrt_2)*sign_R

    # get N
    inv_sqrt_2 = 1./sqrt_2
    dfdy = sign_R*y*(R_T + sqrt_1)/sqrt_1*inv_sqrt_2
    dfdx = sign_R*x*inv_sqrt_2
    inv_absval = 1./np.sqrt(dfdx**2 + dfdy**2 + 1)
    m = inv_absval # multiplier for all values
    if returnformat=='seqrt':
        N = np.column_stack((-dfdx*m, -dfdy*m, m))
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
    else:
        N = (dfdx*m, dfdy*m, m)
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
        #N = (m, dfdx*m, dfdy*m)
        #N = rotate_vector_x(N, surf_rotation=surf_rotation)
    return z, N

def get_z_toric(x, y, R, r, surf_rotation=0, coord_input='original'):
    if coord_input == 'original' and not surf_rotation == 0:
        # rotate the coordinates to align with the cylidner axis
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    R_neg = R < 0
    r_neg = r < 0
    sign_R = 1 - 2*R_neg
    sign_r = 1 - 2*r_neg
    sign_prod = sign_R*sign_r
    # sign_both = 1 - 2*(R_neg and r_neg)
    r = abs(r)
    R = abs(R)
    R_T = R - r*sign_prod

    # get z
    sqrt_1 = sign_prod*np.sqrt(r**2 - y**2)
    sqrt_2 = np.sqrt((R_T + sqrt_1)**2 - x**2)
    return (R - sqrt_2)*sign_R

def get_N_toric(x, y, R, r, surf_rotation=0, coord_input='original', returnformat=''):
    if coord_input == 'original' and not surf_rotation == 0:
        # rotate the coordinates to align with the cylidner axis
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    R_neg = R < 0
    r_neg = r < 0
    sign_R = 1 - 2*R_neg
    sign_r = 1 - 2*r_neg
    sign_prod = sign_R*sign_r
    # sign_both = 1 - 2*(R_neg and r_neg)
    r = abs(r)
    R = abs(R)
    R_T = R - r*sign_prod
    
    sqrt_1 = sign_prod*np.sqrt(r**2 - y**2)
    sqrt_2 = np.sqrt((R_T + sqrt_1)**2 - x**2)

    inv_sqrt_2 = 1./sqrt_2
    dfdy = sign_R*y*(R_T + sqrt_1)/sqrt_1*inv_sqrt_2
    dfdx = sign_R*x*inv_sqrt_2
    inv_absval = 1./np.sqrt(dfdx**2 + dfdy**2 + 1)
    m = inv_absval # multiplier for all values
    if returnformat=='seqrt':
        N = np.column_stack((-dfdx*m, -dfdy*m, m))
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
    else:
        N = (dfdx*m, dfdy*m, m)
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
        #N = (m, dfdx*m, dfdy*m)
        #N = rotate_vector_x(N, surf_rotation=surf_rotation)
    return N