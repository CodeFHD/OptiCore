import numpy as np

from mathutils import Vector

def add_flat_surface(lrad,N1,N2,nsurf=1,xadd=0,yadd=0,nVerts=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """
    
    surfadd=0
    if nsurf == -1:
        surfadd = nVerts+N2-1
    
    verts = []
    faces = []
    
    for j in range(N2)[::nsurf]:
        b = 2*np.pi*j/N2
        verts.append(Vector((-1.*xadd,lrad*np.sin(b)+yadd,lrad*np.cos(b))))
    faces.append([int(surfadd+nsurf*x) for x in range(N2)[::-1]])
    
    return verts, faces