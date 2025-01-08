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

def surftype_zmx2ltype(rx, ry, kx, ky, Ax, Ay):
    """Determine the ltype parameter for the add_lens() function (Blender geometry import)"""
    surf_subtype1 = get_surface_subtype(rx, kx, Ax)
    surf_subtype2 = get_surface_subtype(ry, ky, Ay)
    # evalaute some bools for better readability below
    s1None = surf_subtype1 == 'surf_subtype_None'
    s2None = surf_subtype2 == 'surf_subtype_None'
    s1flat = surf_subtype1 == 'flat'
    s2flat = surf_subtype2 == 'flat'

    # unset None for equality checking
    if rx is None: rx = 0
    if ry is None: ry = 0
    if kx is None: kx = 0
    if ky is None: ky = 0
    if np.all(np.array(Ax) == None): Ax = [0]
    if np.all(np.array(Ay) == None): Ay = [0]
    allequal = rx == ry and kx == ky and Ax == Ay

    # case : bad input
    if s1None and s2None:
        return 'surftype_None'

    # case : flat (return rotational as default)
    elif (s1None and s2flat) or (s2None and s1flat) or (s1flat and s2flat):
        return 'rotational'

    # case : rotational
    elif (not s1None and s2None) or (s1None and not s2None) or allequal:
        return 'rotational'

    # case : cylindricX
    elif not s1None and not s1flat and s2flat:
        return 'cylindricX'

    # case : cylindricY
    elif not s2None and not s2flat and s1flat:
        return 'cylindricY'

    # TODO: not implemented cases (return rotational as default)
    else:
        return 'rotational'


def surftype_Lens2Element(ltype, surf_subtype):
    if surf_subtype == 'flat': return 'flat'
    key = '_'.join((ltype, surf_subtype))
    LOOKUP = {'rotational_spherical': 'spherical',
              'rotational_conic': 'conical',
              'rotational_polynomic': 'polynomic',
              'rotational_aspheric': 'aspheric',
              'cylindricX_spherical': 'cylindrical',
              'cylindricX_conic': 'conicylindric',
              'cylindricX_polynomic': 'polycylindric',
              'cylindricX_aspheric': 'acylindric',
              'cylindricY_spherical': 'cylindrical',
              'cylindricY_conic': 'conicylindric',
              'cylindricY_polynomic': 'polycylindric',
              'cylindricY_aspheric': 'acylindric',}
    """ # not yet implemented
              '', 'toric',
              '', 'conitoric',
              '', 'polytoric',
              '', 'atoric',}
    """
    return LOOKUP[key]

def get_surface_subtype(r, k, A):
    """determine the type of surface w.r.t. the intersection algorithms:"""
    if r is None and np.all(np.array(A) == None):
        return 'surf_subtype_None'

    hasrad = r != 0 and r is not None
    hasconic = k != 0 and k is not None
    haspoly = not (np.all(np.array(A) == 0) or np.all(np.array(A) == None) or A is None)

    if not hasrad and not haspoly:
        surf_subtype = 'flat'
    elif not hasrad and haspoly:
        surf_subtype = 'polynomic'
    elif not hasconic and not haspoly:
        surf_subtype = 'spherical'
    #elif hasrad and not haspoly and k == -1:
    #    surf_subtype = 'parabolic'
    elif hasconic and not haspoly:
        surf_subtype = 'conic'
    else: 
        surf_subtype = 'aspheric'

    return surf_subtype


def check_surface(r, lrad, flrad, k, A, surftype, squarelens):
    ssig = 1 - 2*(r < 0) # this ensures ssig == 1 for r == 0. Sign functions typically return 0.

    surf_subtype = get_surface_subtype(r, k, A)

    if k is None: k = 0

    # get parameters depending on surface type
    if surf_subtype == 'flat':
        # no flange needed for a flat surface
        hasfl = 0
        lrad_surf = lrad
    elif (surftype == 'rotational' and not squarelens) or (surftype != 'rotational' and squarelens):
        hasfl = flrad > 0.01*lrad
        if flrad > 0.99*lrad: flrad = 0.99*lrad
        lrad_surf = lrad - flrad if hasfl else lrad
        if k <= -1:
            # in this case 1+k is smaller than 1 and the situation is safe anyways
            pass
        elif lrad_surf > abs(r)/np.sqrt(1+k):
            lrad_surf = 0.9999*abs(r)/np.sqrt(1+k)
            flrad = lrad - lrad_surf
            hasfl = 1
    else:
        # Flange for other combinations would have non-trivial outline.
        # TODO: This would need significantly more implementation
        # Questioable if this is ever practical
        hasfl = 0
        lrad_surf = lrad
            
    return lrad_surf, hasfl, flrad, ssig, surf_subtype