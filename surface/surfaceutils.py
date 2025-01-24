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
from ..utils.geometry import get_rotmat_z

def get_N1_sqsurface(N2, dshape):
    # passing N2 and calculating N1, not the other way, for historical reasons
    if dshape:
        N1 = N2//2 + 1
    else:
        N1 = N2
    return N1

def get_xy_rotatedsurface(x, y, surf_rotation):
    """
    This function takes x and y coordiantes as input and rotates them around the z-axis,
    one application being rotated surfaces such as cylinder lenses.
    """
    input_is_array = isinstance(x, np.ndarray)
    # rotate the input coordiantes by -surf_rotation to bring them into the local coordinate system
    # R_z_pos = get_rotmat_z(surf_rotation)
    R_z_neg = get_rotmat_z(-surf_rotation)

    if input_is_array:
        arraydim = len(x.shape) 
        if arraydim > 1:
            scp = x.shape
            x = x.ravel()
            y = y.ravel()
        coords_in = np.column_stack((x, y, np.zeros(x.shape[0])))
        coords_rot = np.matmul(R_z_neg, coords_in.T).T
        x = coords_rot[:,0]
        y = coords_rot[:,1]
        if arraydim > 1:
            x = x.reshape(scp)
            y = y.reshape(scp)
    else:
        coords_in = np.array([x, y, 0])
        coords_rot = np.matmul(R_z_neg, coords_in)
        x = coords_rot[0]
        y = coords_rot[1]

    return x, y