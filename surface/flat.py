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

def add_flat_surface(lrad, N1, N2, zadd=0, xadd=0, nVerts=0, hole=False, hrad=0, dshape=False):
    """Flat surface with circular cross-section"""
    
    verts = []
    faces = []

    minb = 0
    maxb = 2*np.pi
    if dshape:
        minb = -np.pi/2
        maxb = np.pi*N2/(N2-1)

    
    if hole:
        for j in range(N2):
            b = maxb*j/N2 + minb
            verts.append([hrad*np.cos(b) + xadd, hrad*np.sin(b), -zadd])
    for j in range(N2):
        b = maxb*j/N2 + minb
        verts.append([lrad*np.cos(b) + xadd, lrad*np.sin(b), -zadd])

    if hole:
        for j in range(N2):
            fi1 = nVerts + (j+1)%N2
            fi2 = nVerts + j
            fi4 = fi1 + N2
            fi3 = fi2 + N2
            faces.append([fi1,fi2,fi3,fi4])
    else:
        faces.append([int(nVerts + x) for x in range(N2)])

    #define normals
    normals = (N2*(1+hole))*[[0, 0, 1]]

    return verts, faces, normals


def add_sqflat_surface(lwidth, N1, N2, zadd=0, xadd=0, nVerts=0, dshape=False):
    """Flat surface with square cross-section"""

    verts = []
    vo = []
    N_tot = 2*(N1+N2) - 4 #number of points around the outline

    if dshape:
        x0 = 0
    else:
        x0 = -lwidth

    """Vertices added in same order as for sqspherical, but outline only"""
    # left side
    for i in range(N2):
        y = 2*(i/(N2-1) - 0.5)*lwidth
        verts.append([x0, y, -zadd])
    # bottom and top, alternating
    for i in range(1, N1-1):
        if dshape:
            x = (i/(N1-1))*lwidth
        else:
            x = 2*(i/(N1-1) - 0.5)*lwidth
        verts.append([x, -lwidth, -zadd])
        verts.append([x, lwidth, -zadd])
    # right side
    for i in range(N2):
        y = 2*(i/(N2-1) - 0.5)*lwidth
        verts.append([lwidth, y, -zadd])

    # Face construction: left-b2t + top-l2r + right-t2b + bottom-r2l
    faces = [i for i in range(N2)] + [2*i + 1 + N2 for i in range(N1-2)] + [i + N2 + 2*(N1 - 2) for i in range(N2)[::-1]] + [2*i + N2 for i in range(N1-2)[::-1]]
    # Face construction: bottom-l2r + right-b2t + top-r2l + left-t2b
    faces = [2*i + N2 for i in range(N1-2)] + [i + N2 + 2*(N1 - 2) for i in range(N2)] + [2*i + 1 + N2 for i in range(N1-2)[::-1]] + [i for i in range(N2)][::-1]

    faces = [[x + nVerts for x in faces]]

    #define normals
    normals = N_tot*[[0, 0, 1]]
                
    # reorder outline verts to go around
    if not dshape:
        vo_row1 = verts[:N1]
        vo_row2 = verts[-N1:]
        vo_col1 = verts[N1:-N1][::2]
        vo_col2 = verts[N1:-N1][1::2]
        vo = vo_row1[::-1] + vo_col1 + vo_row2+ vo_col2[::-1]

    return verts, faces, normals, vo