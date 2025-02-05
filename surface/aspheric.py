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

from .radialprofiles import get_z_evenasphere, get_N_evenasphere, get_dzdr_evenasphere

def _check_k(k, r, lrad):
    #check if k is too large, else crop
    comp = (r/lrad)**2 - 1
    if k > comp*0.9999:
        k = comp*0.9999
    return k



def add_aspheric_surface(R, k, A, lrad, N1, N2, zadd=0, nVerts=0, dshape=False, lrad_ext=0):
    """
    zadd has to be set for second surface (only)
    """

    """
    return add_sagsurface_circular(rad, k, A, lrad, N1, N2,
                                   surftype="aspheric",
                                   zadd=zadd, nVerts=nVerts, dshape=dshape, lrad_ext=lrad_ext)
    """
    
    verts = []
    faces = []
    normals = []

    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N2/(N2-1)

    k = _check_k(k, R, lrad)

    verts.append([-zadd,0,0])
    normals.append((1,0,0))
    r = lrad/(N1- (lrad_ext > lrad))
    x = get_z_evenasphere(r**2, R, k, A)
    for j in range(N2):
        b = maxb*j/N2
        verts.append([-1.*x-zadd, r*np.sin(b), r*np.cos(b)])
        N = get_N_evenasphere(r, b, R, k, A)
        normals.append(N)
        if dshape and j==N2-1:
            pass
        else:
            fi1 = nVerts
            fi2 = fi1 + ((j+1)%N2 + 1)
            fi3 = fi1 + (j + 1)
            faces.append([fi1, fi2, fi3])
    for i in range(1,N1 - (lrad_ext > lrad)):
        r = lrad*(i+1)/(N1 - (lrad_ext > lrad))
        x = get_z_evenasphere(r**2, R, k, A)
        for j in range(N2):
            b = maxb*j/N2
            verts.append([-1.*x-zadd,r*np.sin(b),r*np.cos(b)])
            N = get_N_evenasphere(r, b, R, k, A)
            normals.append(N)
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts + ((j+1) + i*N2)
                fi2 = nVerts + ((j+1)%N2 + 1 + i*N2)
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                faces.append([fi4,fi3,fi2,fi1])
                
    #if there is flat annulus
    if lrad_ext > lrad:
        i = N1 - 2
        r = lrad_ext
        for j in range(N2):
            b = maxb*j/N2
            verts.append([-1.*x - zadd,r*np.sin(b),r*np.cos(b)])
            normals.append((1,0,0))
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts + (j+i*N2+1)+N2
                fi2 = nVerts + ((j+1)%N2+i*N2+1) + N2
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                faces.append([fi4,fi3,fi2,fi1])
            
    return verts, faces, normals

def add_sqaspheric_surface(R, k, A, lwidth, N1, N2, nsurf=1, zadd=0, nVerts=0, cylindrical=False):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    zadd has to be set for second surface (only)
    """

    verts = []
    faces = []
    vertquads = []
    normals = []

    if cylindrical:
        testrad = lwidth/2
    else:
        testrad = lwidth/2*np.sqrt(2)
    k = _check_k(k, R, testrad)

    for i in range(N1):
        y = lwidth*(i/(N1-1) - 0.5)
        for j in range(N2):
            z = lwidth*(j/(N2-1) - 0.5)
            if cylindrical:
                r = y
            else:
                r = np.sqrt(y**2 + z**2)
            x = get_z_evenasphere(r**2, R, k, A)
            verts.append([-1.*x*nsurf-zadd,y,z])
            dxdr = get_dzdr_evenasphere(r, R, k, A)
            adxdr = np.sqrt(dxdr**2 + 1)
            dxdr = dxdr/adxdr
            if cylindrical:
                b = np.pi/2
            else:
                b = np.arctan2(y,z)
            normals.append((nsurf/adxdr, dxdr*np.sin(b), dxdr*np.cos(b)))
            ang = np.arctan2(z+lwidth/(2*N2-2),y+lwidth/(2*N1-2))
            cond1 = N1%2 == 0
            cond2 = j==int(N2/2 - 1) or i==int(N1/2 - 1)
            if cond1 and cond2:
                vertquads.append(-99)
            else:
                vertquads.append(ang/np.pi*2)

    for i in range(N1-1):
        for j in range(N2-1):
            f1 = nVerts + j + 1 + N2*i
            f2 = nVerts + j + N2*i
            f3 = nVerts + j + N2*(i+1)
            f4 = nVerts + j + 1 + N2*(i+1)
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