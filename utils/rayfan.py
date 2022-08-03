import bpy
import numpy as np

def retasnparray(func):
    def wrapper(*args, **kwargs):
        O, D = func(*args, **kwargs)
        O = np.array(O)
        D = np.array(D)
        return O, D
    return wrapper

@retasnparray
def rayfan2D(Nrays, rad, rayfanx=-20):
    #init ray fan (2D plot)
    O = [[rayfanx,0,i] for i in np.linspace(-rad, rad, Nrays)]
    
    D = [[1.,0,0] for i in range(len(O))]
    
    return O, D

@retasnparray
def rayfan3D(Nrays, rad, rayfanx=-20):
    #init ray fan(3D image)
    ao = 0#np.pi/Nrays#angular offset to reduce effect of rays on axis exactly
    O = [[rayfanx, i*np.sin(j+ao), i*np.cos(j+ao)] for i in np.linspace(0,rad,Nrays) for j in np.linspace(0,2*np.pi,Nrays,endpoint=False)]
    #remove duplicate center rays
    O = O[Nrays-1:]
    
    D = [[1.,0,0] for i in range(len(O))]
    
    return O, D

@retasnparray
def rayfan3D_uniformdiskrandom(Nrays, rad, rayfanx=-20):
    #init ray fan (3D uniform sample disk)
    u1 = np.random.rand(Nrays)
    u2 = np.random.rand(Nrays)
    
    Ox = np.zeros(Nrays) + rayfanx
    Oy = rad * np.sqrt(u1) * np.cos(2*np.pi*u2)
    Oz = rad * np.sqrt(u1) * np.sin(2*np.pi*u2)
    
    O = np.vstack((Ox, Oy, Oz)).T
    
    D = [[1.,0,0] for i in range(len(O))]
    
    return O, D

@retasnparray
def rayfan3D_tri(Nrays,rad, rayfanx=-20):
    O = []
    for i in range(Nrays):
        y = rad*(2*i/(Nrays-1) - 1)
        for j in range(Nrays - i%2):
            z = rad*(2*j/(Nrays-1) - 1) + i%2*rad/(Nrays-1)
            O.append([rayfanx,y,z])
    
    D = [[1.,0,0] for i in range(len(O))]
    
    return O, D