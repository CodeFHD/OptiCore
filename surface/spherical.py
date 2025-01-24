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

"""
This _fixN1, _fixN2 mechanic seems to not have been needed
rest of code would have broken due to implementation if it had been used.
Commented out for now, see what happens for some time then delete.

def _fixN1(N1, alpha):
    # recursive function to test and reduce N1 until condition is met
    if (1 - alpha)**(N1 - 2) > 0.01:
        return N1
    else:
        return _fixN1(N1-1, alpha)

def _fixN2(N2, maxb):
    # recursive function to test and increase N2 until condition is met
    if maxb/N2 < 1:
        return N2
    else:
        return _fixN2(N2+1, maxb)
"""



def add_spherical_surface(rad, lrad, N1, N2, zadd=0, nVerts=0, cylinderaxis=None,
                          hole=False, hrad=0, dshape=False, lrad_ext=0):
    #TODO: Could parse arguments by kwargs-dict to keep function call short, then pop here

    verts = []
    faces = []
    normals = []

    minb = 0
    maxb = 2*np.pi
    if dshape:
        minb = -np.pi/2
        maxb = np.pi*N2/(N2-1)

    sig = 1
    if rad < 0:
        sig = -1
    nhole = not hole
    rad = np.abs(rad)
    ang = np.arcsin(lrad/rad)

    # central vertex only without hole
    if not hole:
        verts.append(Vector((0, 0, -zadd)))
        normals.append((0, 0, 1))
        hang = 0 # angle where first ring is placed
    else:
        hang = np.arcsin(hrad/rad)

    # create rings
    for i in range(N1 - (lrad_ext > lrad)):
        a = hang + (ang-hang)*(i+nhole)/(N1 - hole - (lrad_ext > lrad))
        r0 = rad*np.sin(a)
        for j in range(N2):
            b = maxb*j/N2 + minb
            x = r0*np.cos(b)
            y = r0*np.sin(b)
            if cylinderaxis == 'X':
                r = x
            elif cylinderaxis == 'Y':
                r = y
            else:
                r = r0
            z = rad-np.sqrt(rad**2-r**2)
            verts.append(Vector((x, y, -z*sig-zadd)))
            normals.append((sig*np.sin(a)*np.cos(b),
                            sig*np.sin(a)*np.sin(b),
                            np.cos(a)))
            if dshape and j==N2-1:
                pass
            elif i == 0:
                fi1 = nVerts
                fi2 = fi1 + ((j+1)%N2+1)
                fi3 = fi1 + (j+1)
                faces.append([fi3,fi2,fi1])
            else:
                fi1 = nVerts + j+nhole+i*N2
                fi2 = nVerts + (j+1)%N2+nhole+i*N2
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                if cylinderaxis is not None:
                    # in this case faces must be added as tris because 4 vertices will not lie in a plane in general
                    faces.append([fi2, fi3, fi4])
                    faces.append([fi1, fi2, fi4])
                else:
                    faces.append([fi1,fi2,fi3,fi4])
                    #faces.append([fi4,fi3,fi2,fi1])

    #if there is flat annulus, add the outer ring
    if lrad_ext > lrad:
        i = N1 - 2
        r = lrad_ext       
        for j in range(N2):
            b = maxb*j/N2 + minb
            verts.append(Vector((r*np.cos(b), r*np.sin(b), -z*sig - zadd)))
            normals.append((0, 0, 1))
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts + (j+nhole+i*N2)+N2
                fi2 = nVerts + ((j+1)%N2+nhole+i*N2) + N2
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                faces.append([fi1,fi2,fi3,fi4])
                #faces.append([fi4,fi3,fi2,fi1])
            
    return verts, faces, normals#, N1, N2


def add_sqspherical_surface(rad, lwidth, N1, N2, zadd=0, nVerts=0,
                            cylinderaxis=None, dshape=False, lwidth_ext=0):
    """
    zadd has to be set for second surface (only)
    """

    verts = []
    faces = []
    normals = []
    vertquads = []

    sig = 1
    if rad < 0:
        sig = -1
    rad = np.abs(rad)

    hasflangeY = (lwidth_ext > lwidth) and cylinderaxis == 'X'
    hasflangeZ = (lwidth_ext > lwidth) and cylinderaxis == 'Y'

    if cylinderaxis== 'X':
        lwidthY = lwidth
        lwidthZ = lwidth_ext
    elif cylinderaxis== 'Y':
        lwidthY = lwidth_ext
        lwidthZ = lwidth
    else:
        lwidthY = lwidth
        lwidthZ = lwidth

    for i in range(N1):
        if dshape:
            if i==N1-1 and hasflangeY:
                y = lwidth_ext
            else:
                y = (i/(N1 - 1))*lwidthY
        else:
            if i==0 and hasflangeY:
                y = -lwidth_ext
            elif i==N1-1 and hasflangeY:
                y = lwidth_ext
            else:
                y = 2*(i/(N1 - 1) - 0.5)*lwidthY
        for j in range(N2):
            is_flangevert = False
            if j==0 and hasflangeZ:
                z = -lwidthZ
            elif j==N2-1 and hasflangeZ:
                z = lwidthZ
            else:
                z = 2*(j/(N2 - 1) - 0.5)*lwidthZ
            if cylinderaxis == 'X':
                if i==0 and hasflangeY and not dshape:
                    r = 2*((i+1)/(N1 - 1) - 0.5)*lwidthY
                    is_flangevert = True
                elif i==N1-1 and hasflangeY:
                    r = 2*((i-1)/(N1 - 1) - 0.5)*lwidthY
                    is_flangevert = True
                else:
                    r = y
            elif cylinderaxis == 'Y':
                if j==0 and hasflangeZ:
                    r = 2*((j+1)/(N2 - 1) - 0.5)*lwidthZ
                    is_flangevert = True
                elif j==N2-1 and hasflangeZ:
                    r = 2*((j-1)/(N2 - 1) - 0.5)*lwidthZ
                    is_flangevert = True
                else:
                    r = z
            else:
                r = np.sqrt(y**2 + z**2)
            x = rad - np.sqrt(rad**2 - r**2)
            verts.append(Vector((y, z, -x*sig - zadd)))
            ang = np.arctan2(z + lwidth_ext/(N2-1), y + lwidth_ext/(N1-1))
            # Normals
            if is_flangevert:
                normals.append((0, 0, 1))
            else:
                a = np.arcsin(r/rad)
                if cylinderaxis == 'X':
                    b = np.pi/2
                elif cylinderaxis == 'Y':
                    b = 0
                else:
                    b = np.arctan2(y, z)
                normals.append((sig*np.sin(a)*np.sin(b),
                                sig*np.sin(a)*np.cos(b),
                                np.cos(a)))
            # marker for face generation
            cond1 = N2%2 == 0
            cond2 = j==int(N2/2 - 1) or (i==int(N1/2 - 1) and not dshape)
            if cond1 and cond2:
                vertquads.append(-99)
            else:
                vertquads.append(ang/np.pi*2)

    # add the faces of the spherical portion
    for i in range(N1-1):
        for j in range(N2-1):
            f1 = nVerts + j + 1 + N2*i
            f2 = nVerts + j + N2*i
            f3 = nVerts + j + N2*(i+1)
            f4 = nVerts + j + 1 + N2*(i+1)
            if vertquads[j+i*N2] >= 1:
                faces.append([f1,f2,f4])
                faces.append([f2,f3,f4])
            elif vertquads[j+i*N2] >= 0:
                faces.append([f1,f2,f3])
                faces.append([f1,f3,f4])
            elif vertquads[j+i*N2] >= -1:
                faces.append([f1,f2,f4])
                faces.append([f2,f3,f4])
            elif vertquads[j+i*N2] >= -2:
                faces.append([f1,f2,f3])
                faces.append([f1,f3,f4])
            else:
                faces.append([f1,f2,f3,f4])
            
    return verts, faces, normals