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

import os
import re
import numpy as np

from ...raytrace.element import Element
from ...raytrace.lenssystem import Lenssystem
from .surfaces import parse_zmx_surface, get_stop

from ...utils import debugprint

# Because of FLoat rounding, values are rounded to these digits for comparison
EPSILON_DIGITS = 3

def _split_opened_zmx(f):
    lines = f.readlines()
    headerlines = []
    surflines = {} # dict, for each surface index, a list of the lines
    footerlines = []
    lc = 0 # linecounter

    # loop over header
    for line in lines:
        line = line.strip()
        if line.startswith('SURF'):
            break
        headerlines.append(line)
        lc += 1

    # loop over surfaces
    surfidx = -1
    for line in lines[lc:]:
        if re.match(r'\w', line) and not line.startswith('SURF'):
            # break when something else than a surf comes up
            break
        line = line.strip()
        lc += 1
        if line.startswith('SURF'):
            if surfidx > -1:
                surflines[surfidx] = curlines 
            surfidx = int(line.split()[1])
            curlines = []
        else:
            curlines.append(line)
    # finish last surface
    surflines[surfidx] = curlines 

    # loop over footer
    for line in lines[lc:]:
        line = line.strip()
        footerlines.append(line)

    return headerlines, surflines, footerlines

def split_zmx_file(filename):
    # check encoding of file
    f = open(filename, 'rb')
    data = f.read(2)
    if data == b'\xff\xfe':
        encoding = 'utf-16-le'
    else:
        encoding = 'utf-8'
    f.close()

    # split he(ader), su(rfaces), fo(oter)
    with open (filename, 'r', encoding=encoding) as f:    
        he, su, fo = _split_opened_zmx(f)

    return he, su, fo

def identify_elements(surf_infos):
    """
    This function evaluates the previously parsed surfaces and splits it into physical Elements,
    i.e. cemented pieces of glass, or mirrors.
    - Dummy surfaces are discarded (except their CT is counted).
    - Double passes thorugh one Element are detected.
    - The order in which the so created filtered surfaces are passed is determined.
    """
    idx_elements = [] # one list for each element, with the corresponding surface indices
    idx_order = [] # lists the order that surfaces were passed, in case of reflections and multiple passes
    cur_element = []
    direction = True # direction is flipped at mirrors
    
    # some varaible for tracking duplicate surfaces due to reflections, multiple passes, or split-surface-tables
    abspos = 0 # position of the current surface in absolute units
    cur_idx = 0 # count of unique surfaces
    abspos_at_idx = [] # for the unique surfaces: position along the optical axis
    radius_at_idx = [] # for the unique surfaces: radius of curvature
    glass_at_idx = [] # for the unique surfaces: if there is a glass or mirror at that index
        
    # for detection of multi-surface splits, keep track of the last proper surface radius
    abspos_last = -999
    absrad_last = 0
            
    CT_cumulative = {} # cumulative center thickness including dummies. 
    idx_lastvalid = -1
    
    for i in range(1, len(surf_infos)-1): # loop this way just to make sure surfaces are ordered and skip surfaces 0 and -1 (object and image)
        
        debugprint()
        debugprint('IDENTIFY', i)
        debugprint('-'*40)

        glass_prev = surf_infos[i-1]['hasglass']
        glass_here = surf_infos[i+0]['hasglass']
        glass_next = surf_infos[i+1]['hasglass']
        mirror_prev = surf_infos[i-1]['ismirror']
        mirror_here = surf_infos[i+0]['ismirror']
        mirror_next = surf_infos[i+1]['ismirror']
        glasstype_prev = surf_infos[i-1]['glass']
        glasstype_here = surf_infos[i+0]['glass']
        glasstype_next = surf_infos[i+1]['glass']
        
        debugprint('GLASSES', glass_prev, glass_here, glass_next)
        debugprint('GLASSTYPES', glasstype_prev, glasstype_here, glasstype_next)
        debugprint('MIRRORS', mirror_prev, mirror_here, mirror_next)
        
        # TODO: I should probably compare more than just one radius, e.g. in case of toric lens
        # create some compare_surfaces function or something like that to test for duplicates
        absrad_prev = round(abs(surf_infos[i-1]['radius']), EPSILON_DIGITS)
        absrad_here = round(abs(surf_infos[i+0]['radius']), EPSILON_DIGITS)
        absrad_next = round(abs(surf_infos[i+1]['radius']), EPSILON_DIGITS)
        
        abspos_here = round(abspos, EPSILON_DIGITS) # abspos needs to be incremented whatever happens, but the value before increment counts
        abspos = abspos + surf_infos[i]['CT'] # carry without rounding
        abspos_r = round(abspos, EPSILON_DIGITS)
                          
        isdummy = ((glasstype_here == glasstype_prev) or (not glasstype_here and not glass_prev and cur_element == [])) and not mirror_here
        # case 1: The surface is a dummy
        # no glass, nor refraction; no mirror, no reflection
        if isdummy:
            debugprint('DUMMY')
            debugprint(abspos_here, absrad_here, abspos_last, absrad_last, glass_here)
            abspos_last = abspos_here
            if idx_lastvalid != -1:
                CT_cumulative[idx_lastvalid] = CT_cumulative[idx_lastvalid] + surf_infos[i]['CT'] 
            # nothing to do for dummy surfaces
            continue
        # when the dummy check is passed, create the new CT entry
        idx_lastvalid = i
        CT_cumulative[i] = surf_infos[i]['CT'] 
        
        # isopticalsurf = glass_here or glass_prev or mirror_here # condition that the surface is optical
        lastopticalsurf = (glass_prev and not glass_here and not mirror_prev) or mirror_here # condition that the surface ends a stack of optical surfaces. TODO: This assumes that mirrors always end the element, half-mirrors (e.g. polarized pancake lenses) could break this.
        # Check if the following surface would be a duplicate to this one
        split_next = (abspos_r == abspos_here) and (absrad_next == absrad_here)# and glass_here
        split_continued = (abspos_here == abspos_last) and (absrad_here == absrad_last) and i>1
        # Duplicate surface is defined by vertex position along optical axis and absolute radius of curvature
        # TODO, could keep track of reflections and consider sign of radius
        duplicate_here = (abspos_here in abspos_at_idx) and (absrad_here in radius_at_idx) and not split_continued
        
        debugprint(abspos_here, absrad_here, abspos_last, absrad_last, glass_here)
        
        # case 2: The surface is a duplicate of a previous surface.
        if duplicate_here:
            prev_idx = abspos_at_idx.index(abspos_here)
            debugprint('DUPLICATE', prev_idx)
            # case 2.1: The duplicate is replaced because it features glass
            if not glass_at_idx[prev_idx] and direction*glass_here:
                cur_element[-1] = cur_element[-1] + 1 # this condition should only be invoked at the first pass of this surface pair so this should be safe
                glass_at_idx[prev_idx] = True
            # case 2.2: The previous surface is taken as-is
            else:
                idx_order.append(prev_idx + 1) # the tracer later counts starting at surface 1
            abspos_last = abspos_at_idx[prev_idx]
            absrad_last = radius_at_idx[prev_idx]
        # case 3: The surface is new
        else:
            debugprint('NEW', cur_idx)
            # case 3.1: check if the surface is a continuation of a split. If not, finish previous element and start new
            if split_continued:
                debugprint('CONTINUING SPLIT')
                cur_element[-1] = i
                glass_at_idx[-1] = glass_here
                continue
            cur_element.append(i)
            abspos_at_idx.append(abspos_here)
            radius_at_idx.append(absrad_here)
            glass_at_idx.append(glass_here or mirror_here)
            idx_order.append(cur_idx + 1) # the tracer later counts starting at surface 1
            cur_idx = cur_idx + 1
            abspos_last = abspos_here
            absrad_last = absrad_here
            debugprint('LEN', len(abspos_at_idx))
            # case 3.1: It finishes an Element
            if lastopticalsurf and not split_next:
                debugprint('FINISHING', cur_element)
                idx_elements.append(cur_element)
                cur_element = []    
            # case 3.2: It does not finish an Element
            #else:
            #    pass

        if mirror_here:
            direction = not direction
           
    # If there is an unfinished element for the last surface, termiante it
    if cur_element != []:
        debugprint('FINISHING AT END', cur_element)
        idx_elements.append(cur_element)
        
    return idx_elements, idx_order, CT_cumulative

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

def load_from_zmx(filename, testmode=False, logfile=None):
    pixelpitch = 0.01 # default because that is not stored in zmx files AFAIK
    file_path, _ = os.path.split(filename)

    # split he(ader), su(rfaces), fo(oter)
    he, su, fo = split_zmx_file(filename)

    # parse each surface and extract its parameters
    surf_infos = {}
    for idx, surface in su.items():
        surf_info = parse_zmx_surface(surface)
        # make coating filename absolute paths
        if surf_info['coating'][0] in ['DATA']:
            surf_info['coating'][1] = os.path.join(file_path, surf_info['coating'][1])
        surf_infos[idx] = surf_info

    # identify elements, or groups, i.e. cemented pieces
    idx_elements, idx_order, CT_cumulative = identify_elements(surf_infos)
    debugprint()
    debugprint('idx order', idx_order)
    debugprint()

    # create Element() class objects
    elements, dz_sensor = create_elements(surf_infos, idx_elements, CT_cumulative)

    # Create the Lenssystem() class object
    lens = Lenssystem()
    lens.imported_from = filename
    lens.testmode = testmode
    lens.logfile = logfile
    lens.elements = elements
    lens.surf_sequence = idx_order
    lens.num_optical_surfaces = sum([len(l) for l in idx_elements])
    idx_imagesurf = len(surf_infos)-1
    if surf_infos[idx_imagesurf]['hasglass']:
        lens.imagemat = surf_infos[idx_imagesurf]['glass']
    if surf_infos[idx_imagesurf]['rCA'] is not None:
        imleny = 2*surf_infos[idx_imagesurf]['rCA'] # one dimension always stored in rCA
        # second dimension depends on shape
        if surf_infos[idx_imagesurf]['outline_shape'] == 'square':
            imlenx = 2*surf_infos[idx_imagesurf]['rCA_short']
        else:
            imlenx = imleny
    else:
        imlenx = 10.24 # default
        imleny = 10.24 # default
    # convert to number of pixels
    npixx = int(np.ceil(imlenx/pixelpitch))
    npixy = int(np.ceil(imleny/pixelpitch))
    lens.detector['sizex'] = imlenx
    lens.detector['sizey'] = imleny
    lens.detector['npixx'] = npixx
    lens.detector['npixy'] = npixy
    lens.detector['distance'] = dz_sensor
    lens.detector['pixelpitch'] = pixelpitch # default
    stopidx, stoprad, z_stop = get_stop(surf_infos, idx_elements[0][0])
    lens.apertures[0] = {'idx_surface': None, 'z_ap': z_stop, 'shape': 'CIRCLE', 'radius': stoprad, 'n_blades': 0}

    if testmode:
        lens.detector['distance'] = 5000
        lens.build(0.586)

    return lens