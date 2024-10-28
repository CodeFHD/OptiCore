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

from .. import geometry

def rayfan2D_finite(Nrays, rad, rayfanz=-20, rayfany=0):
    #init ray fan (2D plot)
    O1 = [[0, i, 0] for i in np.linspace(-rad, rad, Nrays)]
    O = [[0, rayfany, rayfanz] for i in np.linspace(-rad, rad, Nrays)]
    O, O1 = np.array(O), np.array(O1)
    D = O1 - O
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T
    return O, D

def rayfan2D(Nrays, rad, rayfanz=-20, angle=0):
    #init ray fan (2D plot)
    O = [[0,i, rayfanz] for i in np.linspace(-rad, rad, Nrays)]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if angle != 0:
        rotmat = geometry.get_rotmat_x(angle)
        O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D

def rayfan3D(Nrays, rad, rayfanz=-20, angle=0):
    #init ray fan(3D image)
    ao = np.pi/Nrays#angular offset to reduce effect of rays on axis exactly
    O = [[i*np.sin(j+ao), i*np.cos(j+ao), rayfanz] for i in np.linspace(0,rad,Nrays) for j in np.linspace(0,2*np.pi,Nrays,endpoint=False)]
    #remove duplicate center rays
    O = O[Nrays-1:]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if angle != 0:
        rotmat = geometry.get_rotmat_x(angle)
        O_0 = np.array([[0, O[i,1], 0] for i in range(O.shape[0])])
        # O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D

def rayfan3D_uniformdiskrandom(Nrays, rad, rayfanz=-20, angle=0):
    #init ray fan (3D uniform sample disk)
    u1 = np.random.rand(Nrays)
    u2 = np.random.rand(Nrays)
    Ox = rad * np.sqrt(u1) * np.sin(2*np.pi*u2)
    Oy = rad * np.sqrt(u1) * np.cos(2*np.pi*u2)
    Oz = np.zeros(Nrays) + rayfanz
    O = np.vstack((Ox, Oy, Oz)).T
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if angle != 0:
        rotmat = geometry.get_rotmat_x(angle)
        O_0 = np.array([[0, O[i,1], 0] for i in range(O.shape[0])])
        # O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D

def rayfan3D_square(Nrays, rad, rayfanz=-20, angle=0):
    Nrays = int(np.sqrt(Nrays))
    O = [[rad*(2*j/(Nrays-1) - 1), rad*(2*i/(Nrays-1) - 1), rayfanz] for j in range(Nrays) for i in range(Nrays)]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if angle != 0:
        rotmat = geometry.get_rotmat_x(angle)
        O_0 = np.array([[0, O[i,1], 0] for i in range(O.shape[0])])
        # O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D

def rayfan3D_tri(Nrays,rad, rayfanz=-20, angle=0):
    Nrays = int(np.sqrt(Nrays))
    O = []
    for i in range(Nrays):
        x = rad*(2*i/(Nrays-1) - 1)
        for j in range(Nrays - i%2):
            y = rad*(2*j/(Nrays-1) - 1) + i%2*rad/(Nrays-1)
            O.append([x,y, rayfanz])
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if angle != 0:
        rotmat = geometry.get_rotmat_x(angle)
        O_0 = np.array([[0, O[i,1], 0] for i in range(O.shape[0])])
        # O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D

def rayfan3D_rings(Nrays,rad, rayfanz=-20, angle=0):
    O = []
    N2 = 8 # fixedthe value so that the first ring has 8 points
    N1 = int((1 + np.sqrt(Nrays))/ 2)
    if N1 < 4: N1 = 4 # minimum number of rays for a 3D-distribution
    # center point
    O.append([0,0,rayfanz])
    for i in range(1, N1):
        r = rad*i/(N1 - 1)
        for j in range(i*N2):
            theta = 2*np.pi*j/(i*N2)
            O.append([r*np.sin(theta),r*np.cos(theta), rayfanz])
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    if angle != 0:
        rotmat = geometry.get_rotmat_x(angle)
        O_0 = np.array([[0, O[i,1], 0] for i in range(O.shape[0])])
        # O_0 = np.array([[0,i, 0] for i in np.linspace(-rad, rad, Nrays)])
        O = np.matmul(rotmat, (O-O_0).T).T + O_0
        D = np.matmul(rotmat, D.T).T
    return O, D

# List of implemented distributions
# See in the functions below, to whiuch they correspond
DISTRIBUTIONS = {'2D': rayfan2D,
                 '2D_uniform': rayfan2D,
                 '2D_random': None,
                 '2D_finite': rayfan2D_finite,
                 '3D': rayfan3D_tri,
                 '3D_tri': rayfan3D_tri,
                 '3D_radial': rayfan3D,
                 '3D_square': rayfan3D_square,
                 '3D_random': rayfan3D_uniformdiskrandom,
                 '3D_rings': rayfan3D_rings,
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
                print(f'ERROR: Ray fan distribution "{self.distribution}" is not implemented!')
                return
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
            self.N_history = {}
            self.special_hits = {} # index: history index after which these come.

    def get_rays(self, onlyvalid=True):
        if onlyvalid:
            return self.O[self.idx_hit], self.D[self.idx_hit]
        else:
            return self.O, self.D

    def update(self, O, D, I, OPD, idx_fail, N=None):
            self.O[self.idx_hit] = O
            self.D[self.idx_hit] = D
            self.I[self.idx_hit] = I
            self.OPD[self.idx_hit] = OPD
            if N is not None:
                self.N[self.idx_hit] = N
            if self.store_history:
                self.idx_history = self.idx_history + 1
                self.O_history[self.idx_history] = np.array(self.O)
                if N is not None:
                    self.N_history[self.idx_history] = np.array(self.N)

            self.O[self.idx_hit][idx_fail] = float('nan')
            self.D[self.idx_hit][idx_fail] = float('nan')
            self.I[self.idx_hit][idx_fail] = float('nan')
            self.OPD[self.idx_hit][idx_fail] = float('nan')
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
            print('Error in Rayfan: Called autofocus() without D_tosensor defined.')
            return None
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
