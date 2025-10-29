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

import os
import sys

import numpy as np
pi2 = 2*np.pi # often used

from ..utils import geometry

"""
NOTES ABOUT GEOMETRY:
- theta is interpreted as a rotation about the y-axis, so that nominal displacement is along x-axis!
"""

def rayfan2D(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    # 2D ray fan of source at infinity
    ct, st = np.cos(-theta), np.sin(-theta)
    cp, sp = np.cos(phi), np.sin(phi)
    ca, sa = np.cos(alpha), np.sin(alpha)
    x0, y0, z0 = O_EP
    D = [[st*cp, st*sp, ct] for _ in range(Nrays)]
    D = np.array(D)
    O = [[x0 + r*ca*ct, y0 + r*sa, z0 - r*ca*st] for r in np.linspace(-rad, rad, Nrays, endpoint=True)]
    O = np.array(O)
    t = (rayfanz - O[:,2])/ct # distance ray fan will start beyond EP location.
    O = O + (t*D.T).T
    return O, D

def rayfan2D_finite(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    rayfany = -rayfanz*np.tan(theta)
    O_0 = np.array([[i*np.cos(alpha), i*np.sin(alpha), 0] for i in np.linspace(-rad, rad, Nrays)])
    O = [[rayfany*np.cos(phi), rayfany*np.sin(phi), rayfanz] for i in np.linspace(-rad, rad, Nrays)]
    O, O_0 = np.array(O), np.array(O_0)
    D = O_0 - O
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T
    return O, D

def rayfan2D_aperture_trace(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    # Specifically for reverse aperture tracing.
    # Variables used in the following way for this function:
    # rad: radial position inside aperture
    # rayfanz: aperture location along optical axis
    # theta: half-angle of the fan (should be passed as pi/2 in general)
    # phi: polar angle inside aperture
    # alpha: Rotation of 2D-Fan around optical axis
    ca, sa = np.cos(alpha), np.sin(alpha)
    cp, sp = np.cos(phi), np.sin(phi)
    O = np.array([[rad*cp, rad*sp, rayfanz] for _ in range(Nrays)])
    D = np.array([[ca*np.sin(i), sa*np.sin(i), np.cos(i)] for i in np.linspace(0, theta, Nrays, endpoint=True)])
    return O, D

def rayfan3D_aperture_trace(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    # 3D version of rayfan2D_aperture_trace.
    # Use 2D for rotationally-symmetric lenses, 3D for e.g. anamorphics
    N2 = 8 # fixed value so that the first ring has 8 points
    N1 = int((1 + np.sqrt(Nrays))/ 2)
    if N1 < 4: N1 = 4 # minimum number of rays for a 3D-distribution
    N_total = 1 + int(N1*(N1 - 1)*N2/2)
    cp, sp = np.cos(phi), np.sin(phi)
    O = np.array([[rad*cp, rad*sp, rayfanz] for _ in range(N_total)])
    D = [[0., 0., 1.]]
    for i in range(1, N1):
        t = theta*i/(N1-1)
        ct, st = np.cos(t), np.sin(t)
        d = [[st*np.cos(p + alpha), st*np.sin(p + alpha), ct] for p in np.linspace(0, np.pi, i*N2, endpoint = True)]
        D = D + d
    D = np.array(D)
    return O, D

def rayfan3D_rings(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    N2 = 8 # fixed value so that the first ring has 8 points
    N1 = int((1 + np.sqrt(Nrays))/2)
    if N1 < 4: N1 = 4 # minimum number of rays for a 3D-distribution
    Nrays = 1 + int(N1*(N1 - 1)*N2/2) # replace with new total value
    ct, st = np.cos(-theta), np.sin(-theta)
    cp, sp = np.cos(phi), np.sin(phi)
    # ca, sa = np.cos(alpha), np.sin(alpha)
    x0, y0, z0 = O_EP
    D = [[st*cp, st*sp, ct] for _ in range(Nrays)]
    D = np.array(D)
    O = [O_EP]
    for i in range(1, N1):
        r = rad*i/(N1 - 1)
        x = [r*np.cos(pi2*j/(i*N2) - alpha)*ct for j in range(i*N2)]
        y = [r*np.sin(pi2*j/(i*N2) - alpha) for j in range(i*N2)]
        z = [-1*r*np.cos(pi2*j/(i*N2) - alpha)*st for j in range(i*N2)]
        O_temp = [[x0 + x[j]*cp - y[j]*sp, y0 + x[j]*sp + y[j]*cp, z0 + z[j]] for j in range(i*N2)]
        O += O_temp
    O = np.array(O)
    t = (rayfanz - O[:,2])/ct # distance ray fan will start beyond EP location.
    O = O + (t*D.T).T
    return O, D

def rayfan3D_rings_finite(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    N2 = 8 # fixed value so that the first ring has 8 points
    N1 = int((1 + np.sqrt(Nrays))/ 2)
    if N1 < 4: N1 = 4 # minimum number of rays for a 3D-distribution
    N_total = 1 + int(N1*(N1 - 1)*N2/2) 
    rayfany = -rayfanz*np.tan(theta)
    O = [[rayfany*np.cos(phi), rayfany*np.sin(phi), rayfanz] for i in range(N_total)] # actual origin, all at the same point
    O_0 = []
    O_0.append([0, 0, 0]) # center point
    for i in range(1, N1):
        r = rad*i/(N1 - 1)
        O_temp = [[r*np.sin(pi2*j/(i*N2) - alpha), r*np.cos(pi2*j/(i*N2) - alpha), 0] for j in range(i*N2)]
        O_0 += O_temp
    O, O_0 = np.array(O), np.array(O_0)
    D = O_0 - O
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T
    return O, D

def get_Nrays_trifan(Nrays, rad):
    if Nrays < 25: Nrays = 25 # Minimum for this distribution
    # separate this function to be called in other parts of the code
    Nrays_y = int(np.sqrt(Nrays)) # numhber of rays along x-axis
    Nrays_y = Nrays_y + 1 - Nrays_y%2 # make sure it is an odd number
    a = 2*rad/(Nrays_y - 1) # length of triangle, sampling along x-axis
    h = np.sqrt(3)/2*a # height of triangle, sampling along y-axis
    Nrays_x = int(2*rad//h) + 1 # number of rays along y-axis
    Nrays_x = Nrays_x + 1 - Nrays_x%2 # make sure it is an odd number
    W_half = h*(Nrays_x - 1)/2
    return Nrays_x, Nrays_y, a, h, W_half

def rayfan3D_tri(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    Nrays_x, Nrays_y,a, h, W_half = get_Nrays_trifan(Nrays, rad)
    Nrays = int(Nrays_x*Nrays_y - (Nrays_x-1)//2) # replace with new total value
    # sample lens with equilateral triangles
    ct, st = np.cos(-theta), np.sin(-theta)
    cp, sp = np.cos(phi), np.sin(phi)
    # ca, sa = np.cos(alpha), np.sin(alpha)
    x0, y0, z0 = O_EP
    D = [[st*cp, st*sp, ct] for _ in range(Nrays)]
    D = np.array(D)
    O = []
    for i in range(Nrays_x):
        x = [(i*h - W_half)*ct for j in range(Nrays_y - i%2)]
        y = [(j*a - rad + i%2*a/2) for j in range(Nrays_y - i%2)]
        z = [-1*(i*h - W_half)*st for j in range(Nrays_y - i%2)]
        O_temp = [[x0 + x[j]*cp - y[j]*sp, y0 + x[j]*sp + y[j]*cp, z0 + z[j]] for j in range(Nrays_y - i%2)]
        O += O_temp
    O = np.array(O)
    t = (rayfanz - O[:,2])/ct # distance ray fan will start beyond EP location.
    O = O + (t*D.T).T
    return O, D

def rayfan3D_tri_finite(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    # sample lens with equilateral triangles
    Nrays_x, Nrays_y, a, h, W_half = get_Nrays_trifan(Nrays, rad)
    O_0 = []
    for i in range(Nrays_x):
        x = i*h - W_half
        O_temp = [[x, j*a - rad + i%2*a/2, 0] for j in range(Nrays_y - i%2)]
        O_0 += O_temp
    O_0 = np.array(O_0)
    rayfany = -rayfanz*np.tan(theta)
    N_total = O_0.shape[0]
    O = np.array([[rayfany*np.cos(phi), rayfany*np.sin(phi), rayfanz] for i in range(N_total)])
    if alpha != 0:
        rotmat = geometry.get_rotmat_z(alpha)
        O_0 = np.matmul(rotmat, (O_0).T).T
    D = O_0 - O
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T
    return O, D

def rayfan3D_uniformdiskrandom(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    # alpha makes no sense in random functions, ignoring here
    seed = int.from_bytes(os.urandom(4), sys.byteorder)
    rng = np.random.default_rng(seed)
    u1 = rng.random(Nrays)
    u2 = rng.random(Nrays)
    ct, st = np.cos(-theta), np.sin(-theta)
    cp, sp = np.cos(phi), np.sin(phi)
    # ca, sa = np.cos(alpha), np.sin(alpha)
    x0, y0, z0 = O_EP
    D = [[st*cp, st*sp, ct] for _ in range(Nrays)]
    D = np.array(D)
    Ox = rad * np.sqrt(u1) * np.sin(2*np.pi*u2)
    Oy = rad * np.sqrt(u1) * np.cos(2*np.pi*u2)
    O = [[x0 + Ox[i]*ct*cp - Oy[i]*sp, y0 + Ox[i]*ct*sp + Oy[i]*cp, z0 - Ox[i]*st] for i in range(Nrays)]
    O = np.array(O)
    t = (rayfanz - O[:,2])/ct # distance ray fan will start beyond EP location.
    O = O + (t*D.T).T
    return O, D

def rayfan3D_uniformdiskrandom_finite(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    seed = int.from_bytes(os.urandom(4), sys.byteorder)
    rng = np.random.default_rng(seed)
    u1 = rng.random(Nrays)
    u2 = rng.random(Nrays)
    Ox = rad * np.sqrt(u1) * np.sin(2*np.pi*u2)
    Oy = rad * np.sqrt(u1) * np.cos(2*np.pi*u2)
    rayfany = -rayfanz*np.tan(theta)
    O = [[rayfany*np.cos(phi), rayfany*np.sin(phi), rayfanz] for i in range(Nrays)] # actual origin, all at the same point
    O = np.array(O)
    O_0 = np.vstack((Ox, Oy, np.zeros(Nrays))).T
    D = O_0 - O
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T
    return O, D

"""
def rayfan3D_square(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    Nrays = int(np.sqrt(Nrays))
    O = [[rad*(2*j/(Nrays-1) - 1), rad*(2*i/(Nrays-1) - 1), rayfanz] for j in range(Nrays) for i in range(Nrays)]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if theta != 0:
        rotmat = geometry.get_rotmat_x(theta)
        O_0 = np.array([[0, O[i,1], 0] for i in range(O.shape[0])])
        # O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D
"""

def rayfan3D_sun(Nrays, rad, rayfanz=-20, theta=0, phi=0, alpha=0, O_EP=[0,0,0]):
    # This function is designed to excplicitly sample a 32-arcmin source at infinity
    # independently sample the aperture (as in rayfan3D_uniformdiskrandom)
    # and the direction in a cone
    seed = int.from_bytes(os.urandom(4), sys.byteorder)
    rng = np.random.default_rng(seed)
    ct, st = np.cos(-theta), np.sin(-theta)
    # cp, sp = np.cos(phi), np.sin(phi)
    # ca, sa = np.cos(alpha), np.sin(alpha)
    # x0, y0, z0 = O_EP

    # sampling the cone
    # Using variables T and P as theta and phi inside the cone
    # to distinguish from theta, phi describing the axis of the cone
    u1 = rng.random(Nrays)
    u2 = rng.random(Nrays)
    # sample cos(Theta) linearly
    cosT_max = 0.9999891691779588 # np.cos(32/2/60*np.pi/180)
    cosT = cosT_max + u1*(1 - cosT_max)
    sinT = np.sqrt(1 - cosT*cosT)
    # sample phi directly
    P = u2*2*np.pi
    # convert to points on unit sphere
    D = np.column_stack((np.cos(P)*sinT, np.sin(P)*sinT, cosT))
    # rotate D by theta_phi
    if theta != 0:
        rotmat_theta = geometry.get_rotmat_y(-theta)
        D = np.matmul(rotmat_theta, D.T).T
    if phi != 0:
        rotmat_phi = geometry.get_rotmat_z(phi)
        D = np.matmul(rotmat_phi, D.T).T
    # normalize D to make sure there are no conversion errors left
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T

    # sampling the aperture at z=0
    u1 = rng.random(Nrays)
    u2 = rng.random(Nrays)
    Ox = rad * np.sqrt(u1) * np.sin(2*np.pi*u2)
    Oy = rad * np.sqrt(u1) * np.cos(2*np.pi*u2)
    Oz = np.zeros(Nrays)
    O = np.vstack((Ox, Oy, Oz)).T
    # rotate O by theta_phi
    if theta != 0:
        O = np.matmul(rotmat_theta, O.T).T
    if phi != 0:
        O = np.matmul(rotmat_phi, O.T).T
    # translate O from origin to O_EP
    O = O + np.array(O_EP)

    # project O from the aperture along negative D
    t = (rayfanz - O[:,2])/ct # distance ray fan will start beyond EP location.
    O = O + (t*D.T).T

    return O, D


# List of implemented distributions
# See in the functions below, to whiuch they correspond
DISTRIBUTIONS = {'2D': rayfan2D,
                 '2D_uniform': rayfan2D,
                 '2D_random': None,
                 '2D_finite': rayfan2D_finite,
                 '2D_aperture_trace': rayfan2D_aperture_trace,
                 '3D_aperture_trace': rayfan3D_aperture_trace,
                 '3D': rayfan3D_tri,
                 '3D_tri': rayfan3D_tri,
                 '3D_tri_finite': rayfan3D_tri_finite,
                 '3D_random': rayfan3D_uniformdiskrandom,
                 '3D_random_finite': rayfan3D_uniformdiskrandom_finite,
                 '3D_random_sun': rayfan3D_sun, 
                 '3D_rings': rayfan3D_rings,
                 '3D_rings_finite': rayfan3D_rings_finite,
                 }


class RayFan():
    def __init__(self, distribution='2D', initparams=None, store_history=False):
        if distribution and initparams:
            self.distribution = distribution # Distribution of rays in the fan
            self.initparams = initparams # initialisation parameters for the ray fan
            self.store_history = store_history # if True, stores the origin positions in a dictionary at each update
            self.idx_history = 0
        self.reset()

    def reset(self, distribution=None, initparams=None, store_history=None):
        if distribution and initparams:
            self.distribution = distribution
            self.initparams = initparams
        if store_history:
            self.store_history = store_history
        self.O = [] # Ray origin points
        self.D = [] # Ray direction vectors
        self.D_tosensor = None # Direction of the rays going towards the detector for later focus finding
        self.I = [] # Ray intensities
        self.OPD = [] # Optical path travelled by each ray
        self.idx_hit = [] # Hit or miss marker
        self.N = [] # for storing normal vectors
        if self.distribution:
            if not self.distribution in DISTRIBUTIONS:
                raise ValueError(f'ERROR: Ray fan distribution "{self.distribution}" is not implemented!')
            O, D = DISTRIBUTIONS[self.distribution](*self.initparams)
            self.n_rays = O.shape[0]
            self.O = O
            self.D = D
            self.I = np.ones(self.n_rays)
            self.OPD = np.zeros(self.n_rays)
            self.idx_hit = np.ones(self.n_rays, dtype=bool)
            self.N = np.zeros((self.n_rays, 3))
        if self.store_history:
            self.O_history = {}
            self.O_history[0] = np.array(O)
            #self.N_history = {}
            self.special_hits = {} # index: history index after which these come.

    def get_rays(self, onlyvalid=True):
        if onlyvalid:
            return self.O[self.idx_hit], self.D[self.idx_hit]
        else:
            return self.O, self.D
        
    def get_I(self, onlyvalid=True):
        if onlyvalid:
            return self.I[self.idx_hit]
        else:
            return self.I

    def update(self, O, D, I, OPD, idx_fail, N=None):
        self.O[self.idx_hit] = O
        self.D[self.idx_hit] = D
        if I is not None:
            self.I[self.idx_hit] = I
        if OPD is not None:
            self.OPD[self.idx_hit] = OPD
        if N is not None:
            self.N[self.idx_hit] = N
        if self.store_history:
            self.idx_history = self.idx_history + 1
            self.O_history[self.idx_history] = np.array(self.O)
            #if N is not None:
            #    self.N_history[self.idx_history] = np.array(self.N)
        
        O = self.O[self.idx_hit]
        O[idx_fail] = float('nan')
        self.O[self.idx_hit] = O
        
        D = self.D[self.idx_hit]
        D[idx_fail] = float('nan')
        self.D[self.idx_hit] = D
        
        I = self.I[self.idx_hit]
        I[idx_fail] = float('nan')
        self.I[self.idx_hit] = I
        
        OPD = self.OPD[self.idx_hit]
        OPD[idx_fail] = float('nan')
        self.OPD[self.idx_hit] = OPD
        
        self.idx_hit = ~np.isnan(self.O[:,0])

    def update_special_hits(self, P, idx_fail):
        if not self.store_history:
            return
        idx_special = str(self.idx_history)
        while True:
            if not idx_special in self.special_hits:
                break
            idx_special = idx_special + '.' # simply add a dot to the index to make it different. extract surface using str.split('.')
        self.special_hits[idx_special] = np.full(self.O.shape, float('nan'))
        # print('U', self.special_hits[idx_special].shape, self.idx_hit.shape, idx_fail.shape, P.shape)
        self.special_hits[idx_special][self.idx_hit] = P
        
    def calc_rms_spotsize(self, P=None):
        # get points in detector plane
        if P is None: P = self.O_history[len(self.O_history)-1]
        # get mean point
        Pmean = np.nanmean(P, axis=0)
        # get differences
        Pdiff = P - Pmean
        r2 = np.einsum('ij,ij->i',Pdiff,Pdiff)
        # calculate rms
        rmsspotsize = np.sqrt(np.nanmean(r2))
        return rmsspotsize
    
    def autofocus(self, EFL=None):
        if self.D_tosensor is None:
            raise Exception('Error in Rayfan: Called autofocus() without D_tosensor defined.')
        if EFL is None:
            # if no EFL is given, take BFL from rays
            EFL = 1
        #offsets = np.array([-EFL/100, -EFL/200, 0, EFL/200, EFL/100]) # offset from nominal focus position
        P = self.O_history[len(self.O_history)-1]
        D = self.D_tosensor
        """
        Iteration 1
        """
        offsets = np.array([EFL*i/8 for i in np.linspace(-1,1,31, endpoint=True)]) # in the first pass, coarse scan ca. 10% of the EFL - including margin hence divide by 8
        rms_list = []
        for offset in offsets:
            t = offset/D[:,2] # project offset along the rays
            P_adj = P + (D.T*t).T # compute adjusted points
            rmsspotsize = self.calc_rms_spotsize(P_adj) # compute rms spot size and append
            rms_list.append(rmsspotsize)
        # calculate the focus estimate as the intersection between two lines.
        # skip one point around the minimum to ensure linear region    
        idx = np.where(rms_list == min(rms_list))[0][0]
        if idx < 3 or idx > 27:
            # minimum too close to the edge of the scan range
            return None
        x11, x12, x21, x22 = offsets[idx-3], offsets[idx-2], offsets[idx+3], offsets[idx+2]
        y11, y12, y21, y22 = rms_list[idx-3], rms_list[idx-2], rms_list[idx+3], rms_list[idx+2]
        m1 = (y12 - y11)/(x12 -x11)
        b1 = y11 - m1*x11
        m2 = (y22 - y21)/(x22 -x21)
        b2 = y21 - m2*x21
        offset_min = (b2 - b1)/(m1 - m2)
        
        """
        Iteration 2
        """
        offsets = np.array([offset_min + EFL*i/100 for i in np.linspace(-1,1,51, endpoint=True)]) # in the second pass, scan more finely within 1% of the first iteration estimate
        rms_list = []
        for offset in offsets:
            t = offset/D[:,2] # project offset along the rays
            P_adj = P + (D.T*t).T # compute adjusted points
            rmsspotsize = self.calc_rms_spotsize(P_adj) # compute rms spot size and append
            rms_list.append(rmsspotsize)
        # lowest rms is best focus  
        idx = np.where(rms_list == min(rms_list))[0][0]
        offset_min = offsets[idx]

        t = offset_min/D[:,2]
        # compute adjusted points
        P_new = P + (D.T*t).T
        
        return offset_min, P_new
