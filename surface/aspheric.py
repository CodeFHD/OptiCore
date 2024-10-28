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

def _check_k(k, rad, lrad):
    #check if k is too large, else crop
    comp = (rad/lrad)**2 - 1
    if k > comp*0.9999:
        k = comp*0.9999
    return k

def getdxdr(rad, r, k, A):
    sqt = np.sqrt(1-(1+k)*(r**2/rad**2))
    t1 = 2*r/rad/(1+sqt)
    t2 = (1+k) * r**3
    t22 = rad**3 * sqt * (1 + sqt)**2
    t3 = sum([2*(i+2)*A[i]*r**(2*(i+2) - 1) for i in range(len(A))])
    return t1 + t2/t22 + t3

def getz(rad, r2, k, A):
    f1 = (r2/rad)/(1+np.sqrt(1-(1+k)*(r2/rad**2))) # conic term
    f2 = sum([A[i]*r2**(i+2) for i in range(len(A))]) # polynomic term
    return f1+f2

def get_N(r, b, rad, k, A):
    dxdr = getdxdr(rad, r, k, A)
    adxdr = np.sqrt(dxdr**2 + 1)
    dxdr = dxdr/adxdr
    N = np.column_stack((-dxdr*np.cos(b), -dxdr*np.sin(b), 1/adxdr))
    return N

def add_aspheric_surface(rad, k, A, lrad, N1, N2, nsurf=1, xadd=0, nVerts=0, dshape=False, lrad_ext=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """

    verts = []
    faces = []
    normals = []
    
    surfadd=0
    if nsurf == -1:
        surfadd = N2*N1

    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N2/(N2-1)

    k = _check_k(k, rad, lrad)

    verts.append(Vector((-xadd,0,0)))
    normals.append((nsurf,0,0))
    r = lrad/N1
    x = getz(rad,r**2,k,A)
    for j in range(N2)[::nsurf]:
        b = maxb*j/N2
        verts.append(Vector((-1.*x*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
        dxdr = getdxdr(rad, r, k, A)
        adxdr = np.sqrt(dxdr**2 + 1)
        dxdr = dxdr/adxdr
        normals.append((nsurf/adxdr, dxdr*np.sin(b), dxdr*np.cos(b)))
        if dshape and j==N2-1:
            pass
        else:
            fi1 = nVerts+surfadd
            fi2 = fi1+nsurf*((j+1)%N2+1)
            fi3 = fi1+nsurf*(j+1)
            faces.append([fi1,fi2,fi3])
    for i in range(1,N1 - (lrad_ext > lrad)):
        r = lrad*(i+1)/(N1 - (lrad_ext > lrad))
        x = getz(rad,r**2,k,A)
        for j in range(N2)[::nsurf]:
            b = maxb*j/N2
            verts.append(Vector((-1.*x*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
            dxdr = getdxdr(rad, r, k, A)
            adxdr = np.sqrt(dxdr**2 + 1)
            dxdr = dxdr/adxdr
            normals.append((nsurf/adxdr, dxdr*np.sin(b), dxdr*np.cos(b)))
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts+surfadd+nsurf*(j+1+i*N2)
                fi2 = nVerts+surfadd+nsurf*((j+1)%N2+1+i*N2)
                fi3 = fi2-nsurf*N2
                fi4 = fi1-nsurf*N2
                faces.append([fi4,fi3,fi2,fi1])
                
    #if there is flat annulus
    if lrad_ext > lrad:
        i = N1 - 2
        r = lrad_ext
        for j in range(N2):
            b = maxb*j/N2
            verts.append(Vector((-1.*x - xadd,r*np.sin(b),r*np.cos(b))))
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

def add_sqaspheric_surface(rad,k,A,lwidth,N1,N2,nsurf=1,xadd=0,nVerts=0,cylindrical=False):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """

    verts = []
    faces = []
    vertquads = []
    normals = []
    
    surfadd=0

    if cylindrical:
        testrad = lwidth/2
    else:
        testrad = lwidth/2*np.sqrt(2)
    k = _check_k(k, rad, testrad)

    for i in range(N1):
        y = lwidth*(i/(N1-1) - 0.5)
        for j in range(N2):
            z = lwidth*(j/(N2-1) - 0.5)
            if cylindrical:
                r = y
            else:
                r = np.sqrt(y**2 + z**2)
            x = getz(rad,r**2,k,A)
            verts.append(Vector((-1.*x*nsurf-xadd,y,z)))
            dxdr = getdxdr(rad, r, k, A)
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