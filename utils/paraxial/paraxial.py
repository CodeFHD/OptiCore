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

def refract(y, u, n0, n1, phi):
    u1 = (n0*u - y*phi) / n1
    return y, u1

def transfer(y, u, t):
    y1 = y + u*t
    return y1, u

def reflect(y, u, r):
    u1 = 2./r*y + u
    return y, u1

def trace_lens(y, u, r_list, t_list, n_list, n_elements):
    """
    n_elements is the number of lens elements, e.g. for a singlet n_elements=1.
    n_list contains the refractive indices including ambient medium befoire and after. length must be n_elements+2 (or longer, see comments in code).
    """
    # trace a pair of surface and thickness
    for i in range(n_elements):
        # surface
        r = r_list[i]
        if r == 0:
            phi = 0
        else:
            phi = (n_list[i+1] - n_list[i]) / r
        y, u = refract(y, u, n_list[i], n_list[i+1], phi)
        # thickness
        y, u = transfer(y, u, t_list[i])
    # trace the final surface
    r = r_list[i+1]
    if r == 0:
        phi = 0
    else:
        phi = (n_list[-1] - n_list[i+1]) / r
    y, u = refract(y, u, n_list[i+1], n_list[-1], phi) # n_list[-1] ensures the ambient medium is always taken, so that an overly long list can be taken and the execution covered by n_elements
    return y, u

def calc_BFL(y, u):
    if u == 0:
        return float('inf')
    # t = -y/u
    t = -y/u
    return t

def calc_EFL(y0, u1):
    if u1 == 0:
        return float('inf')
    # technically the same as calc_BFL, just separtate to make it clear the last u and first y are to be used
    t = -y0/u1
    return t