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

modulepath = os.path.dirname(__file__)
glasscatdir = os.path.join(modulepath, 'glasscatalog_data')

N_FALLBACK = 1.6

"""
Section: Catalog imports
"""

filename = os.path.join(glasscatdir, r'glasscatalog_CDGM.pkl')
with open(filename, 'rb') as inputfile:
    coefficients_CDGM = pickle.load(inputfile)

filename = os.path.join(glasscatdir, r'glasscatalog_Hikari.pkl')
with open(filename, 'rb') as inputfile:
    coefficients_Hikari = pickle.load(inputfile)

filename = os.path.join(glasscatdir, r'glasscatalog_Hoya.pkl')
with open(filename, 'rb') as inputfile:
    coefficients_Hoya = pickle.load(inputfile)

filename = os.path.join(glasscatdir, r'glasscatalog_Ohara.pkl')
with open(filename, 'rb') as inputfile:
    coefficients_Ohara = pickle.load(inputfile)

filename = os.path.join(glasscatdir, r'glasscatalog_Schott.pkl')
with open(filename, 'rb') as inputfile:
    coefficients_Schott = pickle.load(inputfile)

filename = os.path.join(glasscatdir, r'glasscatalog_Sumita.pkl')
with open(filename, 'rb') as inputfile:
    coefficients_Sumita = pickle.load(inputfile)


"""
Section: Special definitions e.g. for legacy materials not included in known databases by manufacturers.
When using the coefficient_special dictionary, entries must contain the formula-name as the first element, then coefficients in appropriate format.
"""

if not 'H-ZLAF55A' in coefficients_CDGM:
    coefficients_CDGM['H-ZLAF55A'] = ['Poly-1-4'] + [3.2808637, -1.9066608e-2, 2.8540621e-2, 1.2586912e-3, -5.2847729e-5, 4.2263573e-6]

coefficients_special = {}
coefficients_special['Schott_F7'] = ['tabulated_n', [0.440, 0.4861, 0.5876, 0.6563, 0.700], [1.647073, 1.6378, 1.625358, 1.620207, 1.617707]]


"""
Section: Map catalogs
"""

CATALOG_MAP = {'CDGM': coefficients_CDGM,
               'Hikari': coefficients_Hikari,
               'Hoya': coefficients_Hoya,
               'Ohara': coefficients_Ohara,
               'Schott': coefficients_Schott,
               'Sumita': coefficients_Sumita,
               'special': coefficients_special,}
# Define a default order in whoch glass catalogs are probed.
# The default dict should include all catalogs.
# This default order is rather arbitrary, suggestions why a specific order is to be preferred are welcome.
# At the time of writing this, there are two different materials named F2 and F5 in Schott and CDGM data. So order is important and/or a separate mechnanism to handle catalog-specification is needed.
CATALOG_ORDER_DEFAULT = ['Schott', 'Ohara', 'Hoya', 'CDGM', 'Hikari', 'Sumita', 'special']

# Iterations are made to check for substitutes, e.g. N-prefixed glasses for Schott
SUBSTITUE_PREFIXES = [None,                     # Try without substitution first
                      'N-', 'P-',               # Schott
                      'S-', 'L-',               # Ohara
                      'E-', 'M-', 'MC-', 'MP-', # Hoya
                      'H-', 'D-',               # CDGM
                      'J-', 'Q-',               # Hikari
                      'K-',                     # Sumita
                      ]

"""
Section: Refractive index formulas
"""

def n_Air(wl):
    # Source of equation: https://www.nature.com/articles/srep46111
    B1, B2, C1, C2 = 0.05792105, 0.00167917, 238.0185, 57.362
    n = 1 + B1 / (C1 - wl**2) + B2 / (C2 - wl**2)
    return n

def nSellmeier(wl, C, coefficientorder='sorted'):
    if not len(C)%2 == 0:
        print(f'ERROR: Number of coefficients supplied to nSellmeier() must be even. Returning default value n={N_FALLBACK}')
        return N_FALLBACK
    n_terms = (len(C))//2
    wl2 = wl**2
    if coefficientorder == 'sorted':
        terms = [C[2*i]*wl2/(wl2 - C[2*i+1]) for i in range(n_terms)]
        n2 = 1 + sum(terms)
        return np.sqrt(n2)
    elif coefficientorder == 'separated':
        terms = [C[i]*wl2/(wl2 - C[i+n_terms]) for i in range(n_terms)]
        n2 = 1 + sum(terms)
        return np.sqrt(n2)
    else:
        print(f'ERROR: Invalid argument for coefficientorder in nSellmeier(): {coefficientorder}. Returning default value n={N_FALLBACK}')
        return N_FALLBACK

def nPolynominal(wl, C, npos=1, nneg=1):
    if not len(C) == (1 + npos + nneg):
        print(f'ERROR: Number of coefficients supplied to nPolynominal() must be equal to 1 + npos + nneg. Returning default value n={N_FALLBACK}')
        return N_FALLBACK
    # This is essentially the same as Equation 3 of refractiveindex.info
    constterm = C[0]
    posterm = [C[i + 1]*wl**(2*(i+1)) for i in range(npos)]
    negterm = [C[i + 1 + npos]*wl**(-2*(i+1)) for i in range(nneg)]
    n2 = constterm + sum(posterm) + sum(negterm)
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

# implementation of formulas defined for the refractive-index.info database
def get_n_by_equation(equation_number, C, wl):
    """
    Equations implemented by refractiveindex.info
    """
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
Section: calaog_specific resolving
"""

def get_n_CDGM(wl, C):
    if C[0] == 'Sellmeier2':
        n = nSellmeier(wl, C[1:], coefficientorder='sorted')
        return n
    elif C[0] == 'Poly-1-4':
        n = nPolynominal(wl, C[1:], npos=1, nneg=4)
        return n
    else:
        return -1

def get_n_Hikari(wl, C):
    n = nPolynominal(wl, C, npos=2, nneg=6)
    return n

def get_n_Hoya(wl, C):
    n = nPolynominal(wl, C, npos=1, nneg=4)
    return n

def get_n_Ohara(wl, C):
    n = nSellmeier(wl, C, coefficientorder='separated')
    return n

def get_n_Schott(wl, C):
    n = nSellmeier(wl, C, coefficientorder='separated')
    return n

def get_n_Sumita(wl, C):
    n = nPolynominal(wl, C, npos=1, nneg=4)
    return n

def get_n_tabulated(wl, C):
    print('WARNING: Tabulated refractive index not yet implemented!')
    return -1

CATALOG_RESOLVE_MAP = {'CDGM': get_n_CDGM,
                       'Hikari': get_n_Hikari,
                       'Hoya': get_n_Hoya,
                       'Ohara': get_n_Ohara,
                       'Schott': get_n_Schott,
                       'Sumita': get_n_Sumita,}

""" Parameter-specific  resolve functions """

def _get_t_from_catalog(glassname, catalog_name, wl, substitute_prefix=None):
    if substitute_prefix is not None:
        glassname = substitute_prefix + glassname
    if 't-' + glassname in CATALOG_MAP[catalog_name]:
        t_data = CATALOG_MAP[catalog_name]['t-' + glassname]
        wl_data = [_[0] for _ in t_data]
        t_data = [_[1] for _ in t_data]
        t = np.interp(wl, wl_data, t_data)
        return t
    else:
        return -1

def _get_n_from_catalog(glassname, catalog_name, wl, substitute_prefix=None):
    if substitute_prefix is not None:
        glassname = substitute_prefix + glassname
    if 'n-' + glassname in CATALOG_MAP[catalog_name]:
        coefficients = CATALOG_MAP[catalog_name]['n-' + glassname]
        n = CATALOG_RESOLVE_MAP[catalog_name](wl, coefficients)
        return n
    else:
        return -1

""" Catalog search iteration """

def _iteration_get_value(glassname, wl=0.5875618, substitute_prefix=None, catalog_order=CATALOG_ORDER_DEFAULT, parameter='n'):
    """
    Options for 'parameter':
    'n': refractive idnex
    't': transmission t10
    """
    if not parameter in ['n', 't']:
        raise ValueError("_iteration_get_value() received invalid parameter:", parameter)
    # Case 1: Check if glassname starts with manufacturer string, i.e. explicitly specified
    glassname_prefix = glassname.split('_')[0]
    if glassname_prefix in CATALOG_MAP:
        glassname = '_'.join(glassname.split('_')[1:])
        if parameter == 'n':
            val = _get_n_from_catalog(glassname, glassname_prefix, wl, substitute_prefix=substitute_prefix)
        elif parameter == 't':
            val = _get_t_from_catalog(glassname, glassname_prefix, wl, substitute_prefix=substitute_prefix)
        if val != -1:
            return val

    # Case 2: Search catalogs in specified order
    for catalog_name in catalog_order:
        if parameter == 'n':
            val = _get_n_from_catalog(glassname, catalog_name, wl, substitute_prefix=substitute_prefix)
        elif parameter == 't':
            val = _get_t_from_catalog(glassname, catalog_name, wl, substitute_prefix=substitute_prefix)
        if val != -1:
            return val

    # Case 3: Material not found in this iteration
    return -1

""" Main resolver """

def get_t(glassname, wl=0.5875618, catalog_order=CATALOG_ORDER_DEFAULT, debug_glassname=False):
    """
    Wavelength is given in micrometers
    """

    # catalogs are forced uppercase, do the same here
    glassname = glassname.upper()

    # Case 1: Special cases outside of catalog handling
    if glassname == 'VACUUM':
        # assume no absorption for vacuum
        return 1.
    elif glassname == 'AIR':
        # assume no absorption for air
        return 1.
    elif glassname.startswith('FIXVALUE_'):
        # This needs a method to specify a value in e.g. a zmx file
        return 1.
    elif glassname.startswith('___BLANK'):
        # This needs a method to specify a value in e.g. a zmx file
        return 1.

    # Case 2: Check manufacturer catalogs
    for substitute_prefix in SUBSTITUE_PREFIXES:
        t = _iteration_get_value(glassname, wl, catalog_order = catalog_order,
                            substitute_prefix = substitute_prefix, parameter = 't')
        if t != -1:
            if substitute_prefix is not None:
                glassname = '_'.join(glassname.split('_')[1:])
                print(f'Glas subsitution from {glassname} to {substitute_prefix + glassname}!')
            return t

    # Case 3: Material not found
    # Do nothing here, this case is handled by get_n


def get_n(glassname, wl=0.5875618, catalog_order=CATALOG_ORDER_DEFAULT, debug_glassname=False):
    """
    Wavelength is given in micrometers
    """

    # catalogs are forced uppercase, do the same here
    glassname = glassname.upper()

    # Case 1: Special cases outside of catalog handling
    if glassname == 'VACUUM':
        return 1./n_Air(wl) # Common practise in optics design to specify refractive index relative to (standard-)air rather than vacuum. TODO: Needs double-check if all manufacturers follow the same logic.
        # return 1. # Vacuum-based definition
    elif glassname == 'AIR':
        # It is common that refrctive indices are specified relative to air instead of vacuum
        # For now, treat n_air as 1
        # More explicit confirmation for different sources/mateirals is needed.
        return 1.
        # return n_Air(wl)
    elif glassname.startswith('FIXVALUE_'):
        # fix refractive index
        return float(glassname.split('_')[1])
    elif glassname.startswith('___BLANK'):
        # from .zmx files
        # not 100 percent sure how to interpret this correctly
        # Items 3 and 4 are obviously index and Abbe number, if d- or e-wavelength (or relative to the wavelength specified in the header) unclear, difference is probably marginal. Choosing "d" here because it is probably historically prevalent.
        # Item 1 can eb either 1 or 2, and item 5 may contain small numbers. No idea what they represent.
        nd = float(glassname.split()[3].replace(',', '.'))
        vd = glassname.split()[4]
        # catch odd cases found in available .zmx examples
        if vd == '?':
            return nd
        vd = float(vd.replace(',', '.'))
        if vd == 0:
            return nd
        n = nAbbe_d(wl, [nd, vd])
        return n

    # Case 2: Check manufacturer catalogs
    for substitute_prefix in SUBSTITUE_PREFIXES:
        n = _iteration_get_value(glassname, wl, catalog_order = catalog_order,
                            substitute_prefix = substitute_prefix, parameter = 'n')
        if n != -1:
            if substitute_prefix is not None:
                glassname = '_'.join(glassname.split('_')[1:])
                print(f'Glas subsitution from {glassname} to {substitute_prefix + glassname}!')
            return n

    # Case 3: Material not found
    if debug_glassname:
        return glassname # for automatic analysis of missing glasses
    else:
        print(f'ERROR: Glass type {glassname} not defined in any list! Defaulting to n={N_FALLBACK}...')
        return N_FALLBACK

    
"""
Section: Tests
"""

def test_catalog_duplicates():
    glasses_checked = {}
    duplicate_count = 0
    for catalog in CATALOG_ORDER_DEFAULT:
        current_coefficients = CATALOG_MAP[catalog]
        for glassname in current_coefficients:
            if glassname in glasses_checked:
                print('Duplicate glass found: ', glassname, glasses_checked[glassname], catalog)
                duplicate_count += 1
            else:
                glasses_checked[glassname] = catalog
    print(f'test_catalog_duplicates() found {duplicate_count} duplicates.')

if __name__ == '__main__':
    test_catalog_duplicates()