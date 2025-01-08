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

def get_ringnormals(N, dshape=False):
    normals = []
    
    if dshape:
        maxb = np.pi*N/(N-1)
    else:
        maxb = 2*np.pi

    for j in range(N):
        b = maxb*j/N
        normals.append((0., np.sin(b), np.cos(b)))

    return normals

def get_sqringnormals(N1, N2, dshape=False):
    normals = []

    """Vertices added in same order as for sqspherical, but outline only"""
    # left side
    for i in range(N2):
        normals.append((0., -1., 0.))
    # bottom and top, alternating
    for i in range(N1 - 2):
        normals.append((0., 0., -1.))
        normals.append((0., 0., 1.))
    # right side
    for i in range(N2):
        normals.append((0., 1., 0.))

    return normals