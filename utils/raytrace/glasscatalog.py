import os
import pickle

import numpy as np

# import pickled glass catalog
modulepath = os.path.dirname(__file__)
infname = os.path.join(modulepath, r'glasscatalog_sellmeier.pkl')
with open(infname, 'rb') as infile:
    sellmeier = pickle.load(infile)
infname = os.path.join(modulepath, r'glasscatalog_cauchy.pkl')
with open(infname, 'rb') as infile:
    cauchy = pickle.load(infile)

# special cases
n_vac = 1.0

"""
Conversion functions
"""

def n_Air(wl):
    # https://www.nature.com/articles/srep46111
    B1, B2, C1, C2 = 0.05792105, 0.00167917, 238.0185, 57.362
    n = 1 + B1 / (C1 - wl**2) + B2 / (C2 - wl**2)
    return n

def nSellmeier(wl, A):
    wl2 = wl**2
    n2 = 1 + A[0]*wl2/(wl2 - A[3]) + A[1]*wl2/(wl2 - A[4]) + A[2]*wl2/(wl2 - A[5])
    return np.sqrt(n2)

def nCauchy(wl, A):
    if len(A) < 3:
        print('ERROR: insufficient number of coefficients for nCauchy(). Defualting to n=1.5')
        if isinstance(wl, np.ndarray):
            return np.full(wl.shape, 1.5)
        else:
            return 1.5
    n2 = A[0] + A[1]*wl**2
    for i in range(2, len(A)):
        n2 = n2 + A[i]/wl**(2*(i-1))
    #n2 = A[0] + A[1] * wl ** 2 + A[2] / wl ** 2 + A[3] / wl ** 4 + A[4] / wl ** 6
    if isinstance(n2, np.ndarray):
        n2[n2 < 0] = 0
    elif n2 < 0:
        return 0
    return np.sqrt(n2)

"""
Resolver
"""

def get_n(gname, wl=0.5):
    """
    Wavelength is given in micrometers
    """
    def _iteration_get_n(gname, wl):
        if gname == 'VACUUM':
            return n_vac
        elif gname == 'AIR':
            return n_Air(wl)
        elif gname in sellmeier.keys():
            return nSellmeier(wl, sellmeier[gname])
        elif gname in cauchy.keys():
            return nCauchy(wl, cauchy[gname])
        elif gname.startswith('FIXVALUE_'):
            return float(gname.split('_')[1])
        else:
            return -1
    n = _iteration_get_n(gname, wl)
    if n != -1:
        return n
    else:
        gname_try = 'N-' + gname
    n = _iteration_get_n(gname_try, wl)
    if n != -1:
        print('Glass substitution from', gname, 'to', gname_try)
        return n
    else:
        pass
    print(f'ERROR: Glass type {gname} not defined in any list! Defaulting to n=1.5...')
    return 1.5