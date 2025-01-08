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

from mathutils import Vector

def get_squarefacenormals(n1, n2):
    normals = []

    for i in range(4):
        y = np.sin(i*np.pi/2)
        z = np.cos(i*np.pi/2)
        dn = [(0,y,z) for j in range(n1+n2+2)]
        normals = normals + dn

    return normals

def add_flat_surface(lrad,N1,N2,xadd=0,yadd=0,nVerts=0,hole=False,hrad=0,dshape=False):
    """Flat surface with circular cross-section"""
    
    verts = []
    faces = []

    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N2/(N2-1)
    
    if hole:
        for j in range(N2):
            b = maxb*j/N2
            verts.append(Vector((-1.*xadd,hrad*np.sin(b)+yadd,hrad*np.cos(b))))
    for j in range(N2):
        b = maxb*j/N2
        verts.append(Vector((-1.*xadd,lrad*np.sin(b)+yadd,lrad*np.cos(b))))

    if hole:
        for j in range(N2):
            fi1 = nVerts + (j+1)%N2
            fi2 = nVerts + j
            fi4 = fi1 + N2
            fi3 = fi2 + N2
            faces.append([fi4,fi3,fi2,fi1])
    else:
        faces.append([int(nVerts + x) for x in range(N2)[::-1]])

    fc = np.array(faces).ravel()

    #define normals
    normals = (N2*(1+hole))*[[1,0,0]]

    return verts, faces, normals


def add_sqflat_surface(lwidth, N1, N2, xadd=0,yadd=0,nVerts=0,dshape=False):
    """Flat surface with square cross-section"""

    verts = []

    N_tot = 2*(N1+N2) - 4 #number of points around the outline

    if dshape:
        y0 = 0
    else:
        y0 = -lwidth

    """Vertices added in same order as for sqspherical, but outline only"""
    # left side
    for i in range(N2):
        z = 2*(i/(N2-1) - 0.5)*lwidth
        verts.append(Vector((-1.*xadd, y0, z)))
    # bottom and top, alternating
    for i in range(1, N1-1):
        if dshape:
            y = (i/(N1-1))*lwidth
        else:
            y = 2*(i/(N1-1) - 0.5)*lwidth
        verts.append(Vector((-1.*xadd, y, -lwidth)))
        verts.append(Vector((-1.*xadd, y, lwidth)))
    # right side
    for i in range(N2):
        z = 2*(i/(N2-1) - 0.5)*lwidth
        verts.append(Vector((-1.*xadd, lwidth, z)))

    """
    if dshape:
        verts.append(Vector((-1.*xadd, 0, -lwidth)))
        verts.append(Vector((-1.*xadd, 0, lwidth)))
    else:
        verts.append(Vector((-1.*xadd, -lwidth, -lwidth)))
        verts.append(Vector((-1.*xadd, -lwidth, lwidth)))
    verts.append(Vector((-1.*xadd, lwidth, lwidth)))
    verts.append(Vector((-1.*xadd, lwidth, -lwidth)))
    """

    # Face construction: left-b2t + top-l2r + right-t2b + bottom-r2l
    faces = [i for i in range(N2)] + [2*i + 1 + N2 for i in range(N1-2)] + [i + N2 + 2*(N1 - 2) for i in range(N2)[::-1]] + [2*i + N2 for i in range(N1-2)[::-1]]
    # Face construction: bottom-l2r + right-b2t + top-r2l + left-t2b
    faces = [2*i + N2 for i in range(N1-2)] + [i + N2 + 2*(N1 - 2) for i in range(N2)] + [2*i + 1 + N2 for i in range(N1-2)[::-1]] + [i for i in range(N2)][::-1]

    faces = [[x + nVerts for x in faces]]

    #define normals
    normals = N_tot*[[1,0,0]]

    return verts, faces, normals

def add_rectflat_surface():
    """Flat surface with rectangular cross-section"""
    pass