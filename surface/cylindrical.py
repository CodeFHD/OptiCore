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

from .radialprofiles import get_zN_spherical, get_z_evenasphere, get_N_evenasphere
from .surfaceutils import get_xy_rotatedsurface
from ..utils.geometry import rotate_vector_x, rotate_vector_z

def get_zN_cylindrical(x, y, R, surf_rotation=0, coord_input='original', returnformat=''):
    if coord_input == 'original' and not surf_rotation == 0:
        # rotate the coordinates to align with the cylidner axis
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    r = np.abs(x)
    phi = np.pi*(x<0)
    z, N = get_zN_spherical(r, phi, R, returnformat='')
    if returnformat=='seqrt':
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
    else:
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
        # N = rotate_vector_x(N, surf_rotation=surf_rotation)
    return z, N

def get_z_cylindrical():
    pass

def get_N_cylindrical():
    pass

def get_zN_acylindrical(x, y, R, k=0, A=[0], surf_rotation=0, coord_input='original', returnformat=''):
    if coord_input == 'original' and not surf_rotation == 0:
        # rotate the coordinates to align with the cylidner axis
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    # the nominal orientation of the cylinder (surf_rotation == 0) is along the x-axis
    # hence the radial distance for the cylinder case is along the x-axis
    r = np.abs(x)
    r2 = x**2
    phi = np.pi*(x<0)
    z = get_z_evenasphere(r2, R, k, A)
    N = get_N_evenasphere(r, phi, R, k, A, returnformat=returnformat)
    if returnformat=='seqrt':
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
    else:
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
        # N = rotate_vector_x(N, surf_rotation=surf_rotation)
    return z, N

def get_z_acylindrical(x, y, R, k=0, A=[0], surf_rotation=0, coord_input='original'):
    if coord_input == 'original' and not surf_rotation == 0:
        # rotate the coordinates to align with the cylidner axis
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    # the nominal orientation of the cylinder (surf_rotation == 0) is along the x-axis
    # hence the radial distance for the cylinder case is along the x-axis
    r2 = x**2
    return get_z_evenasphere(r2, R, k, A)

def get_N_acylindrical(x, y, R, k=0, A=[0], surf_rotation=0, coord_input='original', returnformat=''):
    if coord_input == 'original' and not surf_rotation == 0:
        x, y = get_xy_rotatedsurface(x, y, surf_rotation)
    r = np.abs(x)
    phi = np.pi*(x<0)
    N = get_N_evenasphere(r, phi, R, k, A, returnformat=returnformat)
    if returnformat=='seqrt':
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
    else:
        N = rotate_vector_z(N, surf_rotation=surf_rotation)
        # N = rotate_vector_x(N, surf_rotation=surf_rotation)
    return N