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

import numpy as np

""" flat """

def get_zN_flat(r, phi, returnformat=''):
    z = r*np.cos(phi)
    if returnformat=='seqrt':
        num = r.shape[0]
        zeros = np.zeros(num)
        ones = np.ones(num)
        N = np.column_stack((zeros, zeros, ones))
    else:
        N = (0, 0, 1)
        #N = (1, 0, 0)
    return z, N

def get_z_flat(r, phi):
    z = r*np.cos(phi)
    return z

def get_N_flat(r, phi, returnformat=''):
    if returnformat=='seqrt':
        num = r.shape[0]
        zeros = np.zeros(num)
        ones = np.ones(num)
        N = np.column_stack((zeros, zeros, ones))
    else:
        N = (0, 0, 1)
        #N = (1, 0, 0)
    return N

""" spherical """

def get_zN_spherical(r, phi, R, returnformat=''):
    sign_R = 1 - 2*(R < 0)
    R = np.abs(R)
    sqt = np.sqrt(R**2 - r**2)
    z = (R - sqt)*sign_R
    dzdr = r/sqt*sign_R
    adzdr = np.sqrt(dzdr**2 + 1)
    dzdr = dzdr/adzdr
    if returnformat=='seqrt':
        N = np.column_stack((-dzdr*np.sin(phi), -dzdr*np.cos(phi), 1/adzdr))
    else:
        N = (dzdr*np.cos(phi), dzdr*np.sin(phi), 1./adzdr)
        # N = (1./adzdr, dzdr*np.cos(phi), dzdr*np.sin(phi))
    return z, N

def get_z_spherical(r, R):
    sign_R = 1 - 2*(R < 0)
    R = np.abs(R)
    sqt = np.sqrt(R**2 - r**2)
    z = (R - sqt)*sign_R
    return z

def get_N_spherical(r, phi, R, returnformat=''):
    sign_R = 1 - 2*(R < 0)
    sqt = np.sqrt(R**2 - r**2)
    dzdr = r/sqt*sign_R
    adzdr = np.sqrt(dzdr**2 + 1)
    dzdr = dzdr/adzdr
    if returnformat=='seqrt':
        N = np.column_stack((-dzdr*np.sin(phi), -dzdr*np.cos(phi), 1/adzdr))
    else:
        N = (dzdr*np.cos(phi), dzdr*np.sin(phi), 1./adzdr)
        #N = (1./adzdr, dzdr*np.cos(phi), dzdr*np.sin(phi))
    return N

""" conic """

# implement from asphere when/if needed

""" polynomic """

# implement from asphere when/if needed

""" asphere """

def get_dzdr_evenasphere(r, R, k, A):
    sqt = np.sqrt(1 - (1+k)*r**2/R**2)
    term1 = 2*r/R/(1 + sqt) # conic term - derivative of nominator
    term2 = (1 + k) * r**3 # conic term - derivative of denominator pt1
    term2_div = sqt * R**3 * (1 + sqt)**2 # conic term - derivative of denominator pt2
    term3 = sum([2*(i+2)*A[i]*r**(2*(i+2) - 1) for i in range(len(A))])  # polynominal term
    return term1 + term2/term2_div + term3

# def get_zN_asphere(r, phi, R, k, A, returnformat=''):
#    pass

def get_z_evenasphere(r2, R, k, A):
    term1 = r2/R/(1 + np.sqrt(1 - (1+k)*r2/R**2)) # conic term
    term2 = sum([A[i]*r2**(i+2) for i in range(len(A))]) # polynominal term
    return term1 + term2

def get_N_evenasphere(r, phi, R, k, A, returnformat=''):
    dzdr = get_dzdr_evenasphere(r, R, k, A)
    adzdr = np.sqrt(dzdr**2 + 1)
    dzdr = dzdr/adzdr
    if returnformat=='seqrt':
        N = np.column_stack((-dzdr*np.cos(phi), -dzdr*np.sin(phi), 1/adzdr))
    else:
        N = (dzdr*np.cos(phi), dzdr*np.sin(phi), 1./adzdr)
        # N = (1./adzdr, dzdr*np.sin(phi), dzdr*np.cos(phi))
    return N