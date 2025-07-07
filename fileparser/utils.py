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

from ..raytrace.element import Element
from ..utils import debugprint

def get_encoding(filename):
    # check encoding of file
    f = open(filename, 'rb')
    data = f.read(2)
    if data == b'\xff\xfe':
        encoding = 'utf-16-le'
    else:
        encoding = 'utf-8'
    f.close()

    return encoding

def create_elements(surf_infos, idx_elements, CT_cumulative):
    elements = []
    # firstElement = True
    direction = 1 # at a mirror, the direction is flipped and CT is negative
    for i, idx_list in enumerate(idx_elements):
        ele = Element()
        if len(idx_list) == 1:
            if surf_infos[idx_list[0]]['glass'] == 'MIRROR':
                ele.ismirror = True
        # Fill the element data
        for j, idx in enumerate(idx_list):
            ismirror = surf_infos[idx]['glass'] == 'MIRROR'
            if ismirror:
                coating = ['MIRROR', None, None]
            else:
                coating = surf_infos[idx]['coating'] 

            r1 = surf_infos[idx]['radius']
            r2 = surf_infos[idx]['radius2']
            asph1 = surf_infos[idx]['asph']
            asph2 = surf_infos[idx]['asph2']
            surf_rotation = surf_infos[idx]['surf_rotation']

            surftype = surf_infos[idx]['type']
            lenstype = surf_infos[idx]['ltype']
            ele.outline_shape = surf_infos[idx]['outline_shape'] # Assumption that there is no mix between sufaces for one element

            ele.add_surface(surftype, lenstype=lenstype, radius = r1, radius2 = r2, asph = asph1, asph2 = asph2, rCA = surf_infos[idx]['rCA'], CT = CT_cumulative[idx]*direction, material = surf_infos[idx]['glass'], coating = coating, surf_rotation=surf_rotation)
            
        # get the element offsets
        if i == 0:
            dz = 0
        else:
            CT_list = [surf_infos[j]['CT'] for j in range(idx_elements[i-1][-1], idx_list[0])]
            dz = sum(CT_list)
        #if direction == -1:
        #    ele.reverse()
        ele.direction = direction    
        if ele.ismirror:
            direction = -1*direction
        elements.append([ele, dz])
        debugprint()
        debugprint('FINISHED ELEMENT')
        debugprint()
    # get the sensor distance
    CT_list = [surf_infos[j]['CT'] for j in range(idx_elements[-1][-1], len(surf_infos))]
    dz_sensor = sum(CT_list)
    return elements, dz_sensor

def determine_OC_surftype(radX, kX, AX, radY=None, kY=None, AY=None):
    """
    This function takes surface coefficients
    and determines the corresponding surface shape type for further processing.
    The result does not depend on the order of X- and Y-coefficients, but either input is allowed.
    """
    # first, check that either any radius or A coefficients are defined, else return undefined
    noRadX = radX is None
    noRadY = radY is None
    noConicX = kX is None
    noConicY = kY is None
    noPolyX = np.all(np.array(AX) == None) or AX == [] or AX is None
    noPolyY = np.all(np.array(AY) == None) or AY == [] or AY is None
    XisNone = noRadX and noPolyX
    YisNone = noRadY and noPolyY

    # Initial check that anything was passed at all
    if XisNone and YisNone:
        return 'UNDEFINED'

    # In case Y was supplied and X not, flip around.
    # This will simplify case evaluation below
    # because XisNone does not have to be rechecked
    if XisNone:
        noRadX, noConicX, noPolyX, XisNone, noRadY, noConicY, noPolyY, YisNone = noRadY, noConicY, noPolyY, YisNone, noRadX, noConicX, noPolyX, XisNone
        radX, kX, AX, radY, kY, AY = radY, kY, AY, radX, kX, AX

    hasradX = radX != 0 and radX is not None
    hasradY = radY != 0 and radY is not None
    hasconicX = kX != 0 and kX is not None
    hasconicY = kY != 0 and kY is not None
    haspolyX = not (np.all(np.array(AX) == 0) or np.all(np.array(AX) == None) or AX is None)
    haspolyY = not (np.all(np.array(AY) == 0) or np.all(np.array(AY) == None) or AY is None)

    """
    determine flat, rotational, cylindric, or toric case,
      and then subdivide further depending ona vailable coefficients
    flat: flat definition for one or both directions (radius == 0 and no aspheric coefficients)
    rotational: definition given only for one direction, rest is None
    cylindric: any non-flat definition for one direction, flat for the other
    toric: non-flat definition for both
    TODO: freeform

    return types of the above categories:
    Flat options:            flat,
    rotational options:      spherical,    conical,          polynominal,      aspherical,
    cylindrical options:     cylindrical,  conicylindrical,  polycylindrical,  acylindrical,
    toric options:           toric,        conitoric,        polytoric,        atoric
    """

    # first test for flat
    XisFlat = radX == 0 and not haspolyX
    YisFlat = radY == 0 and not haspolyY
    Surfisflat = (XisFlat and YisFlat) or (XisFlat and YisNone)
    if Surfisflat:
        return 'flat'

    # second test for rotational
    # remember: due to intial check-and-flip, XisNone == False is guaranteed here
    IsRotational = YisNone
    if IsRotational:
        if (hasradX and not hasconicX and not haspolyX):
            return 'spherical'
        elif (hasradX and hasconicX and not haspolyX):
            return 'conical'
        elif (not hasradX and haspolyX):
            return 'polynominal'
        else:
            return 'aspheric'

    # third test for cylindric
    isCylindric = YisFlat or (not YisNone and XisFlat)
    # remember: hasradX == True implies XisFlat == False
    if isCylindric:
        if (hasradX and not hasconicX and not haspolyX) or (hasradY and not hasconicY and not haspolyY):
            return 'cylindrical'
        elif (hasradX and hasconicX and not haspolyX) or (hasradY and hasconicY and not haspolyY):
            return 'conicylindrical'
        elif (not hasradX and haspolyX) or (not hasradY and haspolyY):
            return 'polycylindrical'
        else:
            return 'acylindrical'

    # fourth test for toric
    isToric =  (not XisFlat) and (not YisNone and not YisFlat)
    if isToric:
        if (hasradX and not hasconicX and not haspolyX) and (hasradY and not hasconicY and not haspolyY):
            return 'toric'
        elif (hasradX and hasconicX and not haspolyX) and (hasradY and hasconicY and not haspolyY):
            return 'conitoric'
        elif (not hasradX and haspolyX) and (not hasradY and haspolyY):
            return 'polytoric'
        else:
            # currently this covers mixed cases. Are they relevant though, or would that be freeform anyways?
            return 'atoric'

    # Final case (UNDEFINED) that something odd slipped past the above tests
    print('ERROR in function "determine_OC_surftype". Coefficients do not match any implemented surface shape')
    return 'UNDEFINED'