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
# TODO: glass catalogs seem to specify refrafctive indices relativd to air in at least some cases.
# Check and update if needed.
n_vac = 1.0

"""
Conversion functions
"""

def n_Air(wl):
    # Source of equation: https://www.nature.com/articles/srep46111
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

#modern abbe
def calc_cauchyB(ne, ve):
    F = 0.47999
    C = 0.64385
    B = (ne -1)/ve*(F**2*C**2)/(C**2 - F**2)
    return B

def calc_cauchyA(B, ne):
    e = 0.54607
    return ne - B/e**2

def nAbbe(wl, A):
    ne, ve = A[0], A[1]
    B = calc_cauchyB(ne, ve)
    A = calc_cauchyA(B, ne)
    return A + B/wl**2

#old abbe
def calc_cauchyB_nd(nd, vd):
    F = 0.4861327
    C = 0.6562725
    B = (nd -1)/vd*(F**2*C**2)/(C**2 - F**2)
    return B

def calc_cauchyA_nd(B, nd):
    d = 0.5875618
    return nd - B/d**2

def nAbbe_d(wl, A):
    nd, vd = A[0], A[1]
    B = calc_cauchyB_nd(nd, vd)
    A = calc_cauchyA_nd(B, nd)
    return A + B/wl**2

"""
refractiveindex.info database parsing
"""

def get_n_by_equation(equation_number, C, wl):
    if equation_number == 1:
        n_terms = (len(C) - 1)//2
        n2 = 1 + C[0] + np.sum([C[2*i+1]*wl**2/(wl**2 - C[2*i+2]**2) for i in range(n_terms)])
        return np.sqrt(n2)
    elif equation_number == 2:
        n_terms = (len(C) - 1)//2
        n2 = 1 + C[0] + np.sum([C[2*i+1]*wl**2/(wl**2 - C[2*i+2]) for i in range(n_terms)])
        return np.sqrt(n2)
    elif equation_number == 3:
        n_terms = (len(C) - 1)//2
        n2 = C[0] + np.sum([C[2*i+1]*wl**C[2*i+2] for i in range(n_terms)])
        return np.sqrt(n2)
    elif equation_number == 4:
        # this one is quite specific so I'm sticking to a fixed sum
        n2 = C[0] + C[1]*wl**C[2]/(wl**2 - C[3]**C[4]) + C[5]*wl**C[6]/(wl**2 - C[7]**C[8]) + C[9]*wl**C[10] + C[11]*wl**C[12] + C[13]*wl**C[14] + C[15]*wl**C[16]
        return np.sqrt(n2)
    elif equation_number == 5:
        n_terms = (len(C) - 1)//2
        n = C[0] + np.sum([C[2*i+1]*wl**C[2*i+2] for i in range(n_terms)])
        return n
    elif equation_number == 6:
        n_terms = (len(C) - 1)//2
        n = 1 + C[0] + np.sum([C[2*i+1]/(C[2*i+2] - wl**-2) for i in range(n_terms)])
        return n
    elif equation_number == 7:
        # this one is quite specific so I'm sticking to a fixed sum
        n = C[0] + C[1]/(wl**2 - 0.028) + C[2]*(1/(wl**2 - 0.028))**2 + C[3]*wl**2 + C[4]*wl**4 + C[5]*wl**6
        return n
    elif equation_number == 8:
        # this one is quite specific so I'm sticking to a fixed sum
        x = C[0] + C[1]*wl**2/(wl*+2 - C[2]) + C[3]*wl**2
        n2 = (2*x + 1)/(1 - x)
        return np.sqrt(n2)
    elif equation_number == 9:
        # this one is quite specific so I'm sticking to a fixed sum
        n2 = C[0] + C[1]/(wl**2 - C[2]) + C[3]*(wl-C[4])/((wl - C[4])**2 + C[5])
        return np.sqrt(n2)

def get_glass_refractiveindexinfo(gname, catalog=None):
    # riipath = r''
    # catalogfile_nk = os.path.join(riipath, 'catalog-nk.yml')
    pass

"""
Main Resolver
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
        elif gname.startswith('FIXVALUE_'): # fix refractive index
            return float(gname.split('_')[1])
        elif gname.startswith('___BLANK'): # from .zmx files
            # not 100 percent sure how to interpret this correctly
            # Items 3 and 4 are obviously index and Abbe number, if d- or e-wavelength (or relative to the wavelength specified in the header) unclear, difference is probably marginal. Choosing "d" here because it is probably historically prevalent.
            # Item 1 can eb either 1 or 2, and item 5 may contain small numbers. No idea what they represent.
            nd = float(gname.split()[3].replace(',', '.'))
            vd = gname.split()[4]
            # catch odd cases
            if vd == '?':
                return nd
            vd = float(vd.replace(',', '.'))
            if vd == 0:
                return nd
            n = nAbbe_d(wl, [nd, vd])
            return n
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
        gname_try = 'S-' + gname
    n = _iteration_get_n(gname_try, wl)
    if n != -1:
        print('Glass substitution from', gname, 'to', gname_try)
        return n
    else:
        pass
    print(f'ERROR: Glass type {gname} not defined in any list! Defaulting to n=1.5...')
    # return gname # for automatic analysis of missing glasses
    return 1.5