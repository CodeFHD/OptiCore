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

def _fixN1(N1, alpha):
    """
    recursive function to test and reduce N1 until condition is met
    """
    if (1 - alpha)**(N1 - 2) > 0.01:
        return N1
    else:
        return _fixN1(N1-1, alpha)

def _fixN2(N2, maxb):
    """
    recursive function to test and increase N2 until condition is met
    """
    if maxb/N2 < 1:
        return N2
    else:
        return _fixN2(N2+1, maxb)

def get_ringnormals(N, dshape=False, direction=1):
    """
    direction is a multiplier used for inverting the direction in the mirror-hole etc.
    """
    dtp = np.float64
    normals = []

    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N/(N-1)

    for j in range(N):
        b = maxb*j/N
        normals.append((direction*0.,direction*np.sin(b, dtype=dtp),direction*np.cos(b, dtype=dtp)))

    return normals

def add_spherical_surface(rad,lrad,N1,N2,nsurfe=1,xadd=0,nVerts=0,hole=False,hrad=0,dshape=False,optiverts=False,lrad_ext=0):
    #TODO: Could parse arguments by kwargs-dict to keep function call short, then pop here

    dtp = np.float64#for testing

    verts = []
    faces = []
    normals = []

    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N2/(N2-1)

    sig = 1
    if rad < 0:
        sig = -1
    nhole = not hole
    rad = np.abs(rad, dtype=dtp)
    ang = np.arcsin(lrad/rad, dtype=dtp)
    
    if optiverts:
        #for calc of optimised distribution
        N2 = _fixN2(N2, maxb)
        alpha = maxb/N2
        N1 = _fixN1(N1, alpha)
        #create list of rs
        optirs = [lrad]
        for i in range(N1-1):
            optirs.append((1-alpha)*optirs[-1])
        optirs = optirs[::-1]

    if not hole:
        verts.append(Vector((-xadd,0,0)))
        normals.append((1,0,0))
        if optiverts:
            r = optirs[0]
            a = np.arcsin(r/rad, dtype=dtp)
        else:
            a = ang/N1
            r = rad*np.sin(a, dtype=dtp)
        hang = 0
    else:
        hang = np.arcsin(hrad/rad, dtype=dtp)
        r = hrad
        a = hang

    x = rad-np.sqrt(rad**2-r**2, dtype=dtp)
    for j in range(N2):
        b = maxb*j/N2
        verts.append(Vector((-1.*x*sig-xadd,r*np.sin(b, dtype=dtp),r*np.cos(b, dtype=dtp))))
        normals.append((np.cos(a, dtype=dtp),
                        sig*np.sin(a, dtype=dtp)*np.sin(b, dtype=dtp),
                        sig*np.sin(a, dtype=dtp)*np.cos(b, dtype=dtp)))
        if not hole:
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts
                fi2 = fi1 + ((j+1)%N2+1)
                fi3 = fi1 + (j+1)
                faces.append([fi1,fi2,fi3])
    for i in range(1,N1 - (lrad_ext > lrad)):
        if optiverts:
            r = optirs[i]
            a = np.arcsin(r/rad, dtype=dtp)
        else:
            a = hang + (ang-hang)*(i+nhole)/(N1-hole - (lrad_ext > lrad))
            r = rad*np.sin(a, dtype=dtp)
        x = rad-np.sqrt(rad**2-r**2, dtype=dtp)
        for j in range(N2):
            b = maxb*j/N2
            verts.append(Vector((-1.*x*sig-xadd,r*np.sin(b, dtype=dtp),r*np.cos(b, dtype=dtp))))
            normals.append((np.cos(a, dtype=dtp),
                            sig*np.sin(a, dtype=dtp)*np.sin(b, dtype=dtp),
                            sig*np.sin(a, dtype=dtp)*np.cos(b, dtype=dtp)))
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts + j+nhole+i*N2
                fi2 = nVerts + (j+1)%N2+nhole+i*N2
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                faces.append([fi4,fi3,fi2,fi1])

    #if there is flat annulus
    if lrad_ext > lrad:
        i = N1 - 2
        r = lrad_ext       
        for j in range(N2):
            b = maxb*j/N2
            verts.append(Vector((-1.*x*sig - xadd,r*np.sin(b, dtype=dtp),r*np.cos(b, dtype=dtp))))
            normals.append((1,0,0))
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts + (j+nhole+i*N2)+N2
                fi2 = nVerts + ((j+1)%N2+nhole+i*N2) + N2
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                faces.append([fi4,fi3,fi2,fi1])
            
    return verts, faces, normals, N1, N2


def add_sqspherical_surface(rad,lwidth,N1,N2,nsurf=1,xadd=0,nVerts=0,cylindrical=False):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """

    dtp = np.float64

    verts = []
    faces = []
    vertquads = []
    normals = []
    
    surfadd=0

    sig = 1
    if rad < 0:
        sig = -1
    rad = np.abs(rad)

    for i in range(N1):
        y = lwidth*(i/(N1-1) - 0.5)
        for j in range(N2):
            z = lwidth*(j/(N2-1) - 0.5)
            if cylindrical:
                r = y
            else:
                r = np.sqrt(y**2 + z**2)
            x = rad-np.sqrt(rad**2-r**2)
            verts.append(Vector((-1.*x*sig*nsurf-xadd,y,z)))
            ang = np.arctan2(z+lwidth/(2*N2-2),y+lwidth/(2*N1-2))
            a = np.arcsin(r/rad)
            if cylindrical:
                b = np.pi/2
            else:
                b = np.arctan2(y,z)
            normals.append((nsurf*np.cos(a, dtype=dtp),
                            sig*np.sin(a, dtype=dtp)*np.sin(b, dtype=dtp),
                            sig*np.sin(a, dtype=dtp)*np.cos(b, dtype=dtp)))
            cond1 = N1%2 == 0
            cond2 = j==int(N2/2 - 1) or i==int(N1/2 - 1)
            if cond1 and cond2:
                vertquads.append(-99)
            else:
                vertquads.append(ang/np.pi*2)

    for i in range(N1-1):
        for j in range(N2-1):
            f1 = nVerts+surfadd+j+1 + N2*i
            f2 = nVerts+surfadd+j+ N2*i
            f3 = nVerts+surfadd+j + N2*(i+1)
            f4 = nVerts+surfadd+j+1 + N2*(i+1)
            if vertquads[j+i*N2] >= 1:
                faces.append([f1,f2,f4][::nsurf])
                faces.append([f2,f3,f4][::nsurf])
            elif vertquads[j+i*N2] >= 0:
                faces.append([f1,f2,f3][::nsurf])
                faces.append([f1,f3,f4][::nsurf])
            elif vertquads[j+i*N2] >= -1:
                faces.append([f1,f2,f4][::nsurf])
                faces.append([f2,f3,f4][::nsurf])
            elif vertquads[j+i*N2] >= -2:
                faces.append([f1,f2,f3][::nsurf])
                faces.append([f1,f3,f4][::nsurf])
            else:
                faces.append([f1,f2,f3,f4][::nsurf])
            
    return verts, faces, normals