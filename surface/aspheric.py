import numpy as np

from mathutils import Vector

def getz(rad, r,k,A):
    f1 = (r**2/rad)/(1+np.sqrt(1-(1+k)*(r**2/rad**2)))
    f2 = sum([A[i]*r**(2*(i+2)) for i in range(len(A))])
    return f1+f2

def add_aspheric_surface(rad, k, A, lrad, N1, N2,nsurf=1,xadd=0, nVerts=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """

    verts = []
    faces = []
    splitverts = []
    
    surfadd=0
    if nsurf == -1:
        surfadd = N2*N1

    sig = 1
    if rad < 0:
        sig = -1
    rad = np.abs(rad)
    #ang = np.arcsin(lrad/rad)

    verts.append(Vector((-xadd,0,0)))
    splitverts.append(0)
    r = lrad/N1
    x = getz(rad,r,k,A)
    for j in range(N2)[::nsurf]:
        b = 2*np.pi*j/N2
        verts.append(Vector((-1.*x*sig*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
        splitverts.append(0)
        fi1 = nVerts+surfadd
        fi2 = fi1+nsurf*((j+1)%N2+1)
        fi3 = fi1+nsurf*(j+1)
        faces.append([fi1,fi2,fi3])
    for i in range(1,N1):
        r = lrad/N1*(i+1)
        x = getz(rad,r,k,A)
        for j in range(N2)[::nsurf]:
            b = 2*np.pi*j/N2
            verts.append(Vector((-1.*x*sig*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
            if i == N1-1:
                splitverts.append(1)
            else:
               splitverts.append(0)
            fi1 = nVerts+surfadd+nsurf*(j+1+i*N2)
            fi2 = nVerts+surfadd+nsurf*((j+1)%N2+1+i*N2)
            fi3 = fi2-nsurf*N2
            fi4 = fi1-nsurf*N2
            faces.append([fi4,fi3,fi2,fi1])
            
    return verts, faces, splitverts