import numpy as np

def rayfan2D_finite(Nrays, rad, rayfanz=-20, rayfany=0):
    #init ray fan (2D plot)
    O1 = [[0, i, 0] for i in np.linspace(-rad, rad, Nrays)]
    O = [[0, rayfany, rayfanz] for i in np.linspace(-rad, rad, Nrays)]
    O, O1 = np.array(O), np.array(O1)
    D = O1 - O
    D = D.T / np.sqrt(np.einsum('ij,ij->i', D, D)) # normalize
    D = D.T
    return O, D

def rayfan2D(Nrays, rad, rayfanz=-20):
    #init ray fan (2D plot)
    O = [[0,i, rayfanz] for i in np.linspace(-rad, rad, Nrays)]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    return O, D

def rayfan3D(Nrays, rad, rayfanz=-20):
    #init ray fan(3D image)
    ao = np.pi/Nrays#angular offset to reduce effect of rays on axis exactly
    O = [[i*np.sin(j+ao), i*np.cos(j+ao), rayfanz] for i in np.linspace(0,rad,Nrays) for j in np.linspace(0,2*np.pi,Nrays,endpoint=False)]
    #remove duplicate center rays
    O = O[Nrays-1:]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    O, D = np.array(O), np.array(D)
    return O, D

def rayfan3D_uniformdiskrandom(Nrays, rad, rayfanz=-20):
    #init ray fan (3D uniform sample disk)
    u1 = np.random.rand(Nrays)
    u2 = np.random.rand(Nrays)
    Ox = rad * np.sqrt(u1) * np.sin(2*np.pi*u2)
    Oy = rad * np.sqrt(u1) * np.cos(2*np.pi*u2)
    Oz = np.zeros(Nrays) + rayfanz
    O = np.vstack((Ox, Oy, Oz)).T
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    return O, D

def rayfan3D_square(Nrays, rad, rayfanz=-20):
    Nrays = int(np.sqrt(Nrays))
    O = [[rad*(2*j/(Nrays-1) - 1), rad*(2*i/(Nrays-1) - 1), rayfanz] for j in range(Nrays) for i in range(Nrays)]
    D = [[0,0,1.] for i in range(len(O))]
    O, D = np.array(O), np.array(D)
    return O, D

def rayfan3D_tri(Nrays,rad, rayfanz=-20):
    Nrays = int(np.sqrt(Nrays))
    O = []
    for i in range(Nrays):
        x = rad*(2*i/(Nrays-1) - 1)
        for j in range(Nrays - i%2):
            y = rad*(2*j/(Nrays-1) - 1) + i%2*rad/(Nrays-1)
            O.append([x,y, rayfanz])
    D = [[0,0,1.] for i in range(len(O))]
    return np.array(O), np.array(D)

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
                #print(self.O.shape)
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
        print('U', self.special_hits[idx_special].shape, self.idx_hit.shape, idx_fail.shape, P.shape)
        self.special_hits[idx_special][self.idx_hit] = P
