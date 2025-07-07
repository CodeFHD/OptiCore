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

from ..utils import determine_OC_surftype
from ...utils.check_surface import surftype_zmx2ltype

SUPPORTED_SURFTYPES =['STANDARD', 'EVENASPH', 'BICONICX', 'TOROIDAL']


def _type_from_surflines(surflines):
    surftype = 'STANDARD' # default, as zemax file does not need to contain TYPE line (seen in edmund optics doublets)
    for line in surflines:
        # Get the tr
        if line.startswith('TYPE'):
            surftype = line.split()[1]
            # had some files where XASPHERE was used in addition to EVENASPH
            # It seems that in those cases, EVENASPH only goes up to 8th-order coefficients and XASPHERE also to higher ones. Otehr differences I did not find reference to online.
            # For that reason, treating them as the same here.
            if surftype == 'XASPHERE': surftype = 'EVENASPH'
    return surftype

def parse_zmx_surface(surflines):
    # This function returns a dictionary containing the interpreted parameters
    # converted to the OptiCore interpretation.
    # surf_info is filled i subroutines and combined at the end of this function.
    surf_info = {}

    # defaults if parameters not given for this surface
    surftype = None
    isstop = False
    radiusX = None
    kX = None
    AX = []
    radiusY = None
    kY = None
    AY = []
    CT = None
    hasglass = False
    ismirror = False
    outline_shape = 'circular'
    glass = None
    rCA = None
    rCA_short = None
    lrad = None
    coating = ['FRESNEL_0', None, None]

    # get the surface type first in order to know how to interpret the rest
    surftype = _type_from_surflines(surflines)
    if not surftype in SUPPORTED_SURFTYPES:
        raise ValueError(f'ERROR: Zemax surface type "{surftype}" not implemented in OptiCore. Please open a support ticket on GitHub and include a .zmx file ')

    # search the parameters
    for line in surflines:
        # surftype-independent parameters
        if line.startswith('STOP'):
            isstop = True
        elif line.startswith('CURV'):
            curv = float(line.split()[1])
            radiusX = 0 if curv == 0 else 1./curv
        elif line.startswith('DISZ'):
            CT = float(line.split()[1])
        elif line.startswith('GLAS'):
            hasglass = True
            glass = line.split()[1]
            if glass == '___BLANK':
                glass = ' '.join(line.split()[1:])
            if glass == 'MIRROR':
                ismirror = True
                hasglass = False
        elif line.startswith('DIAM'):
            rCA = float(line.split()[1])
        elif line.startswith('MEMA'):
            lrad = float(line.split()[1])
        elif line.startswith('SQAP'):
            outline_shape = 'square'
            rCA_short = float(line.split()[1])
            rCA = float(line.split()[2])
        elif line.startswith('CONI'):
            kX = float(line.split()[1])
        elif line.startswith('OCCT'):
            coat_type = line.split()[1]
            coat_filename = line.split()[2]
            coat_idx = int(line.split()[3])
            coating = [coat_type, coat_filename, coat_idx]

        # surftype-dependent parameters
        elif line.startswith('PARM'):
            iparm = int(line.split()[1])
            if surftype == 'EVENASPH':
                if iparm == 1:
                    # I have yet to see a file that has a power-2 coefficient
                    # (assuming iparm == 1 is for power 2 coefficient)
                    # The OptiCore-geometry function currently starts at power-4, therefore ignoring this for the moment.
                    # TODO: Add power-2 and refactor rest of code to account for this
                    pass
                else:
                    AX.append(float(line.split()[2]))
            elif surftype == 'BICONICX':
                if iparm == 1:
                    radiusY = float(line.split()[2])
            elif surftype == 'TOROIDAL':
                surf_rotation = float(line.split()[2])
        elif line.startswith('XDAT'):
            iparm = int(line.split()[1])
            if surftype == 'TOROIDAL' and iparm == 2:
                radiusY = float(line.split()[2])
    
    # determine the surface type in the OptiCore scheme
    surftype = determine_OC_surftype(radiusX, kX, AX, radiusY, kY, AY)

    # enter default values into undefined variables
    if AX == []: AX = [None]
    if AY == []: AY = [None]
    if rCA_short is None: rCA_short = rCA

    # determine the lenstype for the Lens() class and match conventions
    ltype = surftype_zmx2ltype(radiusX, radiusY, kX, kY, AX, AY)
    if ltype == 'cylindricX':
        ltype = 'cylindrical'
        surf_rotation = np.pi/2
    elif ltype == 'cylindricY':
        radiusX, kX, AX, radiusY, kY, AY = radiusY, kY, AY, radiusX, kX, AX
        ltype = 'cylindrical'
        surf_rotation = 0
    elif ltype == 'toric':
        surf_rotation = 0#np.pi/2
    else:
        surf_rotation = 0

    # fill the return dictionary
    surf_info['type'] = surftype # 'type' is teh OptiCore surftype, i.e. base shape types w/o direction
    surf_info['ltype'] = ltype # ltype is for the Lens class, may include directions, like cylindricalX
    surf_info['isstop'] = isstop
    surf_info['radius'] = radiusX
    surf_info['asph'] = [kX] + AX
    surf_info['radius2'] = radiusY
    surf_info['asph2'] = [kY] + AY
    surf_info['surf_rotation'] = surf_rotation
    #surf_info['radiusX'] = radiusX
    #surf_info['asphX'] = [kX] + AX
    surf_info['CT'] = CT
    surf_info['hasglass'] = hasglass
    surf_info['ismirror'] = ismirror
    surf_info['glass'] = glass
    surf_info['rCA'] = rCA
    surf_info['rCA_short'] = rCA_short
    surf_info['lrad'] = lrad
    surf_info['outline_shape'] = outline_shape
    surf_info['coating'] = coating
    return surf_info

def get_stop(surf_infos, idx_first=0):
    stopidx = 1
    for idx, surf in surf_infos.items(): # Use first surface as default stop location
        if idx == 1:
            if surf['rCA'] is not None:
                stoprad = surf['rCA']
            elif surf['lrad'] is not None:
                stoprad = surf['lrad']
            break
    for idx, surf in surf_infos.items(): # override stop where 'isstop' is set
        if surf['isstop']:
            stopidx = idx
            if surf['rCA'] is not None:
                stoprad = surf['rCA']
            elif surf['lrad'] is not None:
                stoprad = surf['lrad']
    if idx_first > stopidx: # Aperture in front of first glass surface
        CT_list = [surf_infos[i]['CT'] for i in range(stopidx, idx_first)]
        z_stop = -sum(CT_list)
    else:
        CT_list = [surf_infos[i]['CT'] for i in range(idx_first, stopidx)]
        z_stop = sum(CT_list)
    return stopidx, stoprad, z_stop