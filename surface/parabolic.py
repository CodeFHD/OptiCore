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

def getdzdr(r,A):
    return 2.*r*A

def add_parabolic_surface(fp, mrad, N1, N2, theta, orig='FP', zadd=0, nVerts=0, hole=False, hrad=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    zadd has to be set for second surface (only)

    orig: FP=focal point, MC = mirror center
    """

    verts = []
    faces = []
    normals = []

    fp *= -1
    nhole = not hole

    #compute basic paramters
    ct = np.cos(theta*2*np.pi/360)
    st = np.sin(theta*2*np.pi/360)
    fs = 2*fp/(1 + ct)
    A = 1/(4*fp)
    OAD = fs*st
    zOAD = OAD**2*A
    zmin = 0

    if orig=='FP':
        zoffset = -fp
        xoffset = 0
        fzoffset = -fp + zmin
        fxoffset = OAD
    elif orig == 'MC':
        zoffset = -zOAD
        xoffset = -OAD
        fzoffset = -zOAD + zmin
        fxoffset = 0

    if not hole:
        verts.append([OAD + xoffset, 0, A*OAD**2 + zoffset - zadd])
        dzdr = getdzdr(OAD, A)
        adzdr = np.sqrt(1 + dzdr**2)
        dzdr = dzdr/adzdr
        b = np.pi/2
        normals.append((-1*dzdr*np.sin(b), -1*dzdr*np.cos(b), 1/adzdr))

    for j in range(N2):
        ri = mrad/N1
        if hole:
            ri = hrad
        tj = 2*np.pi*j/N2
        ctj = np.cos(tj)
        stj = np.sin(tj)
        xp = OAD + ri*ctj
        yp = ri*stj
        dp = np.sqrt(xp**2 + yp**2)
        zp = A*dp**2
        verts.append([xp + xoffset, yp, zp + zoffset - zadd])
        dzdr = getdzdr(dp,A)
        adzdr = np.sqrt(1 + dzdr**2)
        dzdr = dzdr/adzdr
        sinb = xp/dp
        cosb = yp/dp
        normals.append((-1*dzdr*sinb, -1*dzdr*cosb, 1/adzdr))
        if not hole:
            fi1 = nVerts
            fi2 = fi1 + (j+1)%N2 + 1
            fi3 = fi1 + j + 1
            faces.append([fi3,fi2,fi1])

    for i in range(1, N1):
        ri = mrad*(i+1)/N1
        if hole:
            ri = hrad + (mrad - hrad)*i/(N1 - 1)
        for j in range(N2):
            tj = 2*np.pi*j/N2
            ctj = np.cos(tj)
            stj = np.sin(tj)
            xp = OAD + ri*ctj
            yp = ri*stj
            dp = np.sqrt(xp**2 + yp**2)
            zp = A*dp**2
            verts.append([xp + xoffset, yp, zp + zoffset - zadd])
            dzdr = getdzdr(dp, A)
            adzdr = np.sqrt(1 + dzdr**2)
            dzdr = dzdr/adzdr
            sinb = xp/dp
            cosb = yp/dp
            normals.append((-1*dzdr*sinb, -1*dzdr*cosb, 1/adzdr))
            fi1 = nVerts + nhole + i*N2 + j
            fi2 = nVerts + nhole + i*N2 + (j+1)%N2
            fi3 = fi2 - N2
            fi4 = fi1 - N2
            faces.append([fi1,fi2,fi3,fi4])

    zs = [v[2] for v in verts]
    fzoffset = min(zs)

    return verts, faces, fxoffset, fzoffset, normals