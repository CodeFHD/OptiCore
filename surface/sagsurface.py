import numpy as np

from mathutils import Vector

from .radialprofiles import get_zN_spherical, get_z_spherical, get_N_spherical, get_z_evenasphere, get_N_evenasphere
# from .spherical import get_z_spherical, get_N_spherical
from .aspheric import _check_k # get_z_asphere, get_N_asphere, 
from .toric import get_zN_toric
from .cylindrical import get_zN_cylindrical, get_zN_acylindrical


def add_sagsurface_circular(R, k, A, lrad, N1, N2,
                            surftype='aspheric',
                            rad2=None, k2=None, A2=None,
                            surf_rotation=0,
                            zadd=0, nVerts=0, dshape=False, lrad_ext=0):
    verts = []
    faces = []
    normals = []

    make_triface = surftype in ['cylindrical', 'acylindrical', 'toric']

    minb = 0
    maxb = 2*np.pi
    if dshape:
        minb = -np.pi/2
        maxb = np.pi*N2/(N2-1)

    if surftype=='aspheric':
        k = _check_k(k, R, lrad)

    verts.append(Vector((0, 0, -zadd)))
    normals.append((0, 0, 1))
    """outer loop"""
    for i in range(N1 - (lrad_ext > lrad)):
        r = lrad*(i+1)/(N1 - (lrad_ext > lrad))
        # for rotational surfaces, get vertex in outer loop for efficiency
        if surftype == 'spherical':
            z = get_z_spherical(r, R)
        elif surftype == 'aspheric':
            z =  get_z_evenasphere(r**2, R, k, A)
        """inner loop"""
        for j in range(N2):
            phi = maxb*j/N2 + minb
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            # for rotational surfaces, get only normal, else get both vertex and normal
            if surftype == 'spherical':
                N = get_N_spherical(r, phi, R)
            elif surftype == 'aspheric':
                N = get_N_evenasphere(r, phi, R, k, A)
            elif surftype == 'cylindrical':
                z, N = get_zN_cylindrical(x, y, R, surf_rotation=surf_rotation)
            elif surftype == 'acylindrical':
                z, N = get_zN_acylindrical(x, y, R, k, A, surf_rotation=surf_rotation)
            elif surftype == 'toric':
                z, N = get_zN_toric(x, y, R, rad2, surf_rotation=surf_rotation)
            verts.append(Vector((x, y, -z - zadd)))
            normals.append(N)
            if dshape and j==N2-1:
                pass
            elif i==0:
                # between center vertex and first ring, there are always triangles
                fi1 = nVerts
                fi2 = fi1+((j+1)%N2+1)
                fi3 = fi1+(j+1)
                faces.append([fi3,fi2,fi1])
            else:
                fi1 = nVerts+(j+1+i*N2)
                fi2 = nVerts+((j+1)%N2+1+i*N2)
                fi3 = fi2-N2
                fi4 = fi1-N2
                if make_triface:
                    # for surfaces without rotational symmetry, need tris
                    faces.append([fi2, fi3, fi4])
                    faces.append([fi1, fi2, fi4])
                else:
                    # for surfaces with rotational symmetry, can make quads
                    faces.append([fi1,fi2,fi3,fi4])
   
    #
    # DO NOT ADD CODE HERE WITHOUT OBSERVING THE NOTE BELOW
    #

    # if there is flat annulus
    # ATTENTION: this uses values from the last loop above.
    # Take care not to overwrite if code is added.
    if lrad_ext > lrad:
        i = N1 - 2
        r = lrad_ext
        for j in range(N2):
            phi = maxb*j/N2 + minb
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            verts.append(Vector((x, y, -z - zadd)))
            normals.append((0, 0, 1))
            if dshape and j==N2-1:
                pass
            else:
                fi1 = nVerts + (j+i*N2+1)+N2
                fi2 = nVerts + ((j+1)%N2+i*N2+1) + N2
                fi3 = fi2 - N2
                fi4 = fi1 - N2
                faces.append([fi1,fi2,fi3,fi4])

    return verts, faces, normals

def add_sagsurface_rectangular(R, k, A, lwidth, N1, N2,
                            surftype='aspheric',
                            rad2=None, k2=None, A2=None,
                            surf_rotation=0,
                            zadd=0, nVerts=0, dshape=False, lrad_ext=0):
    verts = []
    faces = []
    vertquads = []
    normals = []

    if surftype=='aspheric':
        testrad = lwidth*np.sqrt(2)
        k = _check_k(k, R, testrad)

    for i in range(N1):
        if dshape:
            x = (i/(N1 - 1))*lwidth
        else:
            x = 2*(i/(N1 - 1) - 0.5)*lwidth
        for j in range(N2):
            y = 2*(j/(N2 - 1) - 0.5)*lwidth
            r = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            if surftype == 'spherical':
                z, N = get_zN_spherical(r, phi, R)
                # z =  get_z_spherical(r, R)
                # N = get_N_spherical(r, phi, R)
            elif surftype == 'aspheric':
                z = get_z_evenasphere(r**2, R, k, A)
                N = get_N_evenasphere(r, phi, R, k, A)
            elif surftype == 'cylindrical':
                z, N = get_zN_cylindrical(x, y, R, surf_rotation=surf_rotation)
            elif surftype == 'acylindrical':
                z, N = get_zN_acylindrical(x, y, R, k, A, surf_rotation=surf_rotation)
            elif surftype == 'toric':
                z, N = get_zN_toric(x, y, R, rad2, surf_rotation=surf_rotation)
            verts.append(Vector((x, y, -z-zadd)))
            normals.append(N)
            ang = np.arctan2(y + lwidth/(2*N2-2), x + lwidth/(2*N1-2))
            vertquads.append(ang/np.pi*2)

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