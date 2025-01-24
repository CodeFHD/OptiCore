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
        return 'flat'

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
              'rotational_polynominal': 'polynominal',
              'rotational_aspheric': 'aspheric',
              'cylindrical_spherical': 'cylindrical',
              'cylindrical_conic': 'conicylindrical',
              'cylindrical_polynominal': 'polycylindrical',
              'cylindrical_aspheric': 'acylindrical',
              'toric_spherical': 'toric'}
    """ # not yet implemented
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
        surf_subtype = 'flat' # 'polynominal' # pure polynominal not yet supported
    elif not hasconic and not haspoly:
        surf_subtype = 'spherical'
    #elif hasrad and not haspoly and k == -1:
    #    surf_subtype = 'parabolic'
    elif hasconic and not haspoly:
        surf_subtype = 'conic'
    else: 
        surf_subtype = 'aspheric'

    return surf_subtype

def check_surface(lrad, flrad, RX, kX, AX, RY, kY, AY, surftype, squarelens, rd=0):
    # rd is short for recursiondepth
    if rd > 2:
        print("ERROR: check_surface() appears to be caught in a loop!!!")
        return None
    
    # special case of flat surface will be used in multiple cases
    returnvalues_flat = (lrad, 0, 0, 'flat', 'flat', 'flat', 0, 0, [0], 0)

    if kX is None: kX = 0
    if kY is None: kY = 0
    if AX is None or AX == [] or np.all(np.array(AX) == [None]): AX = [0]
    if AY is None or AY == [] or np.all(np.array(AY) == [None]): AY = [0]
    surf_subtype_X = get_surface_subtype(RX, kX, AX)
    surf_subtype_Y = get_surface_subtype(RY, kY, AY)

    # No flange unless explicitly confirmed
    hasfl = 0
    lrad_surf = lrad
    # same for surface rotation offset
    dSurfrot = 0

    """ Case 1: Flat """
    if surftype == 'flat' or (surf_subtype_X == 'flat' and surf_subtype_Y == 'flat'):
        return returnvalues_flat

    """ Case 2: Rotational """
    # all Y values are ignored for this case
    if surftype == 'rotational':
        if surf_subtype_X == 'flat':
            return returnvalues_flat
        if not squarelens:
            if flrad > 0.99*lrad: # maximum flange width
                return returnvalues_flat
            hasfl = flrad > 0.01*lrad # minimum flange width
            lrad_surf = lrad - flrad if hasfl else lrad
            if kX <= -1: # for k<-1, the argument in the sqrt of the conic term is always positive, i.e. well-defined
                pass
            elif lrad_surf > abs(RX)/np.sqrt(1+kX): # automatic flange to maximum of well-defined radius
                lrad_surf = 0.9999*abs(RX)/np.sqrt(1+kX)
                flrad = lrad - lrad_surf
                hasfl = 1

    """ Case 3: Cylindrical """
    # all Y values are ignored for this case
    if surftype == 'cylindrical':
        if surf_subtype_X == 'flat':
            return returnvalues_flat
        if squarelens:
            if flrad > 0.99*lrad: # maximum flange width
                return returnvalues_flat
            hasfl = flrad > 0.01*lrad # minimum flange width
            lrad_surf = lrad - flrad if hasfl else lrad
            if kX <= -1: # for k<-1, the argument in the sqrt of the conic term is always positive, i.e. well-defined
                pass
            elif lrad_surf > abs(RX)/np.sqrt(1+kX): # automatic flange to maximum of well-defined radius
                lrad_surf = 0.9999*abs(RX)/np.sqrt(1+kX)
                flrad = lrad - lrad_surf
                hasfl = 1

    """ Case 4: Toric """
    if surftype == 'toric':
        if surf_subtype_X == 'flat' and surf_subtype_Y == 'flat':
            return returnvalues_flat
        if surf_subtype_Y == 'flat' and not surf_subtype_X == 'flat':
            # equivalent to a cylinder, X-axis already valid
            surftype = 'cylindrical'
            # perform a recursion because this will resolve the proper evaluation of the cylinder surface
            # return check_surface(lrad, flrad, RX, kX, AX, RY, kY, AY, surftype, squarelens, rd=rd+1)
        elif surf_subtype_X == 'flat' and not surf_subtype_Y == 'flat':
            # equivalent to a cylinder, flip parameters
            surftype = 'cylindrical'
            surf_subtype_X, surf_subtype_Y = surf_subtype_Y, surf_subtype_X
            RX, kX, AX, RY, kY, AY = RY, kY, AY, RX, kX, AX
            dSurfrot = np.pi/2
            # return check_surface(lrad, flrad, RX, kX, AX, RY, kY, AY, surftype, squarelens, rd=rd+1)

    """ Default return case """
    return lrad_surf, hasfl, flrad, surftype, surf_subtype_X, surf_subtype_Y, RX, kX, AX, dSurfrot