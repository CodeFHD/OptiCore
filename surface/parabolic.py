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

def getdxdr(r,A):
    return 2.*r*A

def add_parabolic_surface(fp,mrad,N1,N2,theta,orig='FP',nsurf=1,xadd=0,nVerts=0,hole=False,hrad=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)

    orig: FP=focal point, MC = mirror center
    """

    verts = []
    faces = []
    normals = []

    surfadd=0
    if nsurf == -1:
        surfadd = N2*N1

    fp *= -1
    nhole = not hole

    #compute basic paramters
    ct = np.cos(theta*2*np.pi/360)
    st = np.sin(theta*2*np.pi/360)
    fs = 2*fp/(1+ct)
    A = 1/(4*fp)
    OAD = fs*st
    xOAD = OAD**2*A
    xmin = 0

    if orig=='FP':
        xoffset = -fp
        yoffset = 0
        fxoffset = -fp+xmin
        fyoffset = OAD
    elif orig == 'MC':
        xoffset = -xOAD
        yoffset = -OAD
        fxoffset = -xOAD+xmin
        fyoffset = 0

    if not hole:
        verts.append(Vector((A*OAD**2+xoffset-xadd,OAD+yoffset,0)))
        dxdr = getdxdr(OAD,A)
        adxdr = np.sqrt(1 + dxdr**2)
        dxdr = dxdr/adxdr
        b = np.pi/2
        normals.append((1/adxdr, -1*dxdr*np.sin(b), -1*dxdr*np.cos(b)))

    for j in range(N2)[::nsurf]:
        ri = mrad/N1
        if hole:
            ri = hrad
        tj = 2*np.pi*j/N2
        ctj = np.cos(tj)
        stj = np.sin(tj)
        yp = OAD + ri*stj
        zp = ri*ctj
        dp = np.sqrt(yp**2+zp**2)
        xp = A*dp**2
        verts.append(Vector((xp+xoffset-xadd,yp+yoffset,zp)))
        dxdr = getdxdr(dp,A)
        adxdr = np.sqrt(1 + dxdr**2)
        dxdr = dxdr/adxdr
        sinb = yp/dp
        cosb = zp/dp
        normals.append((1/adxdr, -1*dxdr*sinb, -1*dxdr*cosb))
        if not hole:
            fi1 = nVerts+surfadd
            fi2 = fi1+nsurf*((j+1)%N2+1)
            fi3 = fi1+nsurf*(j+1)
            faces.append([fi1,fi2,fi3])

    for i in range(1,N1):
        ri = mrad*(i+1)/N1
        if hole:
            ri = hrad + (mrad-hrad)*(i)/(N1-1)
        for j in range(N2)[::nsurf]:
            tj = 2*np.pi*j/N2
            ctj = np.cos(tj)
            stj = np.sin(tj)
            yp = OAD + ri*stj
            zp = ri*ctj
            dp = np.sqrt(yp**2+zp**2)
            xp = A*dp**2
            verts.append(Vector((xp+xoffset-xadd,yp+yoffset,zp)))
            dxdr = getdxdr(dp,A)
            adxdr = np.sqrt(1 + dxdr**2)
            dxdr = dxdr/adxdr
            sinb = yp/dp
            cosb = zp/dp
            normals.append((1/adxdr, -1*dxdr*sinb, -1*dxdr*cosb))
            fi1 = nVerts+surfadd+nsurf*(j+nhole+i*N2)
            fi2 = nVerts+surfadd+nsurf*((j+1)%N2+nhole+i*N2)
            fi3 = fi2-nsurf*N2
            fi4 = fi1-nsurf*N2
            faces.append([fi4,fi3,fi2,fi1])

    xs = [v.x for v in verts]
    fxoffset = min(xs)

    return verts, faces, fyoffset, fxoffset, normals