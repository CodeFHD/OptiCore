import numpy as np

from mathutils import Vector

def add_spherical_surface(rad,lrad,N1,N2,nsurf=1,xadd=0,nVerts=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """

    verts = []
    faces = []
    
    surfadd=0
    if nsurf == -1:
        surfadd = N2*N1

    sig = 1
    if rad < 0:
        sig = -1
    rad = np.abs(rad)
    ang = np.arcsin(lrad/rad)

    verts.append(Vector((-xadd,0,0)))
    a = ang/N1
    r = rad*np.sin(a)
    x = rad-np.sqrt(rad**2-r**2)
    for j in range(N2)[::nsurf]:
        b = 2*np.pi*j/N2
        verts.append(Vector((-1.*x*sig*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
        fi1 = nVerts+surfadd
        fi2 = fi1+nsurf*((j+1)%N2+1)
        fi3 = fi1+nsurf*(j+1)
        faces.append([fi1,fi2,fi3])
    for i in range(1,N1):
        a = ang/N1*(i+1)
        r = rad*np.sin(a)
        x = rad-np.sqrt(rad**2-r**2)
        for j in range(N2)[::nsurf]:
            b = 2*np.pi*j/N2
            verts.append(Vector((-1.*x*sig*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
            fi1 = nVerts+surfadd+nsurf*(j+1+i*N2)
            fi2 = nVerts+surfadd+nsurf*((j+1)%N2+1+i*N2)
            fi3 = fi2-nsurf*N2
            fi4 = fi1-nsurf*N2
            faces.append([fi4,fi3,fi2,fi1])
            
    return verts, faces