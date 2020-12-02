import numpy as np

from mathutils import Vector

def get_squarefacenormals(n1, n2):
    dtp = np.float64
    normals = []

    for i in range(4):
        y = np.sin(i*np.pi/2)
        z = np.cos(i*np.pi/2)
        dn = [(0,y,z) for j in range(n1+n2+2)]
        normals = normals + dn

    return normals

def add_flat_surface(lrad,N1,N2,nsurf=1,xadd=0,yadd=0,nVerts=0,hole=False,hrad=0,dshape=False):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """
    
    verts = []
    faces = []
    
    surfadd=0
    if nsurf == -1:
        surfadd = N2*(1+hole)-1

    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N2/(N2-1)
    
    if hole:
        for j in range(N2)[::nsurf]:
            b = maxb*j/N2
            verts.append(Vector((-1.*xadd,hrad*np.sin(b)+yadd,hrad*np.cos(b))))
    for j in range(N2)[::nsurf]:
        b = maxb*j/N2
        verts.append(Vector((-1.*xadd,lrad*np.sin(b)+yadd,lrad*np.cos(b))))
    if hole:
        for j in range(N2)[::nsurf]:
            fi1 = nVerts+surfadd+nsurf*((j+1)%N2)
            fi2 = nVerts+surfadd+nsurf*(j)
            fi4 = fi1+nsurf*N2
            fi3 = fi2+nsurf*N2
            faces.append([fi4,fi3,fi2,fi1])
    else:
        faces.append([int(nVerts+surfadd+nsurf*x) for x in range(N2)[::-1]])
    fc = np.array(faces).ravel()

    #define normals
    normals = (N2*(1+hole))*[[nsurf,0,0]]

    return verts, faces, normals

def add_sqflat_surface(lwidth,N1,N2,nsurf=1,xadd=0,yadd=0,nVerts=0):
    
    verts = []
    faces = [[nVerts+1, nVerts, nVerts+3, nVerts+2]]

    verts.append(Vector((-1.*xadd,-0.5*lwidth,-0.5*lwidth)))
    verts.append(Vector((-1.*xadd,-0.5*lwidth,0.5*lwidth)))
    verts.append(Vector((-1.*xadd,0.5*lwidth,0.5*lwidth)))
    verts.append(Vector((-1.*xadd,0.5*lwidth,-0.5*lwidth)))

    #define normals
    normals = 4*[[nsurf,0,0]]

    return verts, faces, normals