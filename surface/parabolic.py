import numpy as np

from mathutils import Vector

def add_parabolic_surface(fp,mrad,N1,N2,theta,orig='FP',nsurf=1,xadd=0,nVerts=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)

    orig: FP=focal point, MC = mirror center
    """

    verts = []
    faces = []

    surfadd=0
    if nsurf == -1:
        surfadd = N2*N1

    fp *= -1

    #compute basic paramters
    ct = np.cos(theta*2*np.pi/360)
    st = np.sin(theta*2*np.pi/360)
    fs = 2*fp/(1+ct)
    A = 1/(4*fp)
    OAD = fs*st
    xOAD = OAD**2*A
    xmin = 0

    if orig=='FP':
        if fp < 1:
            xmin = (np.abs(OAD)+mrad)**2*A
        elif np.abs(OAD) > mrad:
            xmin = (np.abs(OAD)-mrad)**2*A
        xoffset = -fp
        yoffset = 0
        fxoffset = -fp+xmin
        fyoffset = OAD
    elif orig == 'MC':
        if fp < 1:
            xmin = (np.abs(OAD)+mrad)**2*A
        elif np.abs(OAD) > mrad:
            xmin = (np.abs(OAD)-mrad)**2*A
        xoffset = -xOAD
        yoffset = -OAD
        fxoffset = -xOAD+xmin
        fyoffset = 0

    verts.append(Vector((A*OAD**2+xoffset-xadd,OAD+yoffset,0)))

    for j in range(N2)[::nsurf]:
        ri = mrad/N1
        tj = 2*np.pi*j/N2
        ctj = np.cos(tj)
        stj = np.sin(tj)
        yp = OAD + ri*stj
        zp = ri*ctj
        dp = np.sqrt(yp**2+zp**2)
        xp = A*dp**2
        verts.append(Vector((xp+xoffset-xadd,yp+yoffset,zp)))
        fi1 = nVerts+surfadd
        fi2 = fi1+nsurf*((j+1)%N2+1)
        fi3 = fi1+nsurf*(j+1)
        faces.append([fi1,fi2,fi3])

    for i in range(1,N1):
        ri = mrad*(i+1)/N1
        for j in range(N2)[::nsurf]:
            tj = 2*np.pi*j/N2
            ctj = np.cos(tj)
            stj = np.sin(tj)
            yp = OAD + ri*stj
            zp = ri*ctj
            dp = np.sqrt(yp**2+zp**2)
            xp = A*dp**2
            verts.append(Vector((xp+xoffset-xadd,yp+yoffset,zp)))
            fi1 = nVerts+surfadd+nsurf*(j+1+i*N2)
            fi2 = nVerts+surfadd+nsurf*((j+1)%N2+1+i*N2)
            fi3 = fi2-nsurf*N2
            fi4 = fi1-nsurf*N2
            faces.append([fi4,fi3,fi2,fi1])

    return verts, faces, fyoffset, fxoffset