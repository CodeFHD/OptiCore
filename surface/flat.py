import numpy as np

from mathutils import Vector

def add_flat_surface(lrad,N1,N2,nsurf=1,xadd=0,yadd=0,nVerts=0,hole=False,hrad=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """
    
    surfadd=0
    if nsurf == -1:
        surfadd = N2*(1+hole)-1
    
    verts = []
    faces = []
    splitverts = []
    
    if hole:
        for j in range(N2)[::nsurf]:
            b = 2*np.pi*j/N2
            verts.append(Vector((-1.*xadd,hrad*np.sin(b)+yadd,hrad*np.cos(b))))
            splitverts.append(0)
    for j in range(N2)[::nsurf]:
        b = 2*np.pi*j/N2
        verts.append(Vector((-1.*xadd,lrad*np.sin(b)+yadd,lrad*np.cos(b))))
        splitverts.append(1)
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

    return verts, faces, splitverts

def add_sqflat_surface(lwidth,N1,N2,nsurf=1,xadd=0,yadd=0,nVerts=0):
    
    verts = []
    faces = [[nVerts+1, nVerts, nVerts+3, nVerts+2]]
    splitverts = [1,1,1,1]

    verts.append(Vector((-1.*xadd,-0.5*lwidth,-0.5*lwidth)))
    verts.append(Vector((-1.*xadd,-0.5*lwidth,0.5*lwidth)))
    verts.append(Vector((-1.*xadd,0.5*lwidth,0.5*lwidth)))
    verts.append(Vector((-1.*xadd,0.5*lwidth,-0.5*lwidth)))

    return verts, faces, splitverts