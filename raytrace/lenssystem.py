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

from . import glasscatalog
from . import paraxial
from .coatingdata import CoatingData

"""
TODO: Possibly to be foreseen in this implementation:
In case of a multi-bounce system, i.e. where the same surface is used multiple times,
implement a logic that it is not necessary to create the same element twice.
"""


class Lenssystem():
    def __init__(self, name='NONAME', ambimat='AIR'):
        self.name = name
        self.imported_from = None

        # "construction" of the system, i.e. physical contents
        self.elements = [] # list of elements comprising this lens [[Element, displacement]]
        self.apertures = {} # dict of apertures ['idx_surface' (aperture comes after this surface), 'z_ap' (absolute position), 'shape', 'radius', 'n_blades']
        self.detector = {} # detector definition ['distance' from last lens, 'sizex', 'sizey', 'npixx', 'npixy', 'pixelpitch', 'pixelsize']

        # optical parameters
        self.wl = 0.546074 # default wavelength is e-line. unit micrometers
        self.ambimat = ambimat # material between elements and, by default, image and object medium
        self.imagemat = None # optional material behind last lens. Defaults to ambimat.
        self.objectmat = None # optional material before first lens. Defaults to ambimat.
        self.surf_sequence = None # The nominal sequence of surface indices. If None, each surface is passed one, in order. This can be used for multi-bounce systems or ghosts.

        self.idx_image = -1 # index of surface where image is rendered. default -1 (last surface). Used if rendering intermediate images

        # the final data structure
        self.data = {}
        self.coating_data = {}
        self.num_surfaces = 0
        self.optical_axis_current = [0.,0.,1.] # For folded system, this variable keeps track of the local optical axis after the addition of each element. The first element is oriented, by definition, along the z-axis.
        
        # state parameters
        self.isreverse = False
        self.isbuilt = False
        
        # for testing purposes
        self.testmode = False
        self.logfile = None
    
    def _add_coating_data(self, coating):
        # determine key for self.coating_data
        if coating[0] == 'MIRROR':
            key = 'MIRROR'
        elif coating[0] == 'DATA':
            filename = os.path.basename(coating[1])
            key = f'DATA_{filename}'
        elif coating[0] == 'FIXVALUE':
            key = f'FIXVALUE_{coating[1]:.4f}'
        else:
            key = 'FRESNEL_0'
        # create new CoatingData instance if needed
        if not key in self.coating_data:
            new_coating = CoatingData(coating[0])
            if coating[0] == 'DATA':
                new_coating.load_coating_csv(coating[1])
            elif coating[0] == 'FIXVALUE':
                new_coating.reflectance_data = coating[1]
            self.coating_data[key] = new_coating
        return key

    def add_element(self, ele, dz=0, position=None, rotation=[0,0,0], rotation_absolute=False):
        """
        position: defines the [x, y, z] coordinates of the first surface vertex
        rotation: defines the TBD-rotation to describe the element rotation relative to the global coordiante system
        rotation_absolute: determines if the rotation is taken relative to the current, cumulative rotation, or absolute. The optical axis will be 
        if position==None, the element will be added along self.optical_axis_current with the offset dz.
        """
        self.elements.append([ele, dz])

    def build(self, wl=0.5786):
        """This function builds the final data structure"""
        # Data structure definition
        # This dictionary is what the ray tracer actually works thorugh
        self.data['type'] = [] # Surface type for specific functions
        self.data['radius'] = [] # Spherical surface radius [r0, r1, ...]. 
        self.data['asph'] = [] # Apsheric coefficients [[k, a4, a6,...], ...].
        self.data['radius2'] = [] # secondary radius list for toric lenses
        self.data['asph2'] = [] # secondary aspheric list for toric lenses
        # self.data['freeform_coeff'] = [] # coefficients for freeform-equation
        # self.data['CT'] = [] # Center thickness
        self.data['CT_sum'] = [] # Center thickness sum
        self.data['lrad'] = [] # Mechanical Lens radius
        self.data['rCA'] = [] # Clear aperture radius
        self.data['rCA_short'] = [] # Clear aperture radius
        self.data['surf_rotation'] = [] # angle by which surface is rotated around optical axis. E.g. for cylinder lenses.
        self.data['CA_shape'] = [] # Shape of clear aperture: circular, square, or rectangular
        self.data['surf_decenter'] = [] # Surface lateral position error [[x, y]]
        self.data['surf_tilt'] = [] # Surface tilt error [[X, Y]]
        self.data['n'] = [] # Refractive index
        self.data['ismirror'] = [] # If the surface is reflective
        self.data['coating'] = [] # coatings [key, data]
        # self.data['stop'] = [] # Aperture stops [[x,y,z], shape, size]
        # Other parameters
        if wl:
            self.wl = wl
        else:
            self.wl = 0.5876

        surfacecount = 0
        # First element of structure is ambient medium plus dummy surface parameters
        if self.objectmat:
            n = glasscatalog.get_n(self.objectmat, self.wl)
        else:
            n = glasscatalog.get_n(self.ambimat, self.wl)
        self.data['type'].append('OBJECT')
        self.data['radius'].append(None)
        self.data['asph'].append([None, None])
        self.data['radius2'].append(None)
        self.data['asph2'].append([None, None])
        # self.data['freeform_coeff'].append([None])
        # self.data['CT'].append(None)
        self.data['CT_sum'].append(None)
        self.data['lrad'].append(None)
        self.data['rCA'].append(None)
        self.data['rCA_short'].append(None)
        self.data['surf_rotation'].append(None)
        self.data['CA_shape'].append(None)
        self.data['surf_decenter'].append(None)
        self.data['surf_tilt'].append(None)
        self.data['n'].append(n)
        self.data['ismirror'].append(False)
        self.data['coating'].append([None])
        surfacecount = surfacecount + 1
        # Add glass surfaces
        CT_sum = 0
        maxrad = 0
        # lastdz = 0
        for ele, dz in self.elements:
            num_surfaces = len(ele.data['radius'])
            max_rCA = max(ele.data['rCA'])
            for s in range(num_surfaces):
                if s == 0:
                    CT_sum = CT_sum + dz
                else:
                    CT_sum = CT_sum + ele.data['CT'][s-1]
                if s == num_surfaces-1:
                    n = glasscatalog.get_n(self.ambimat, self.wl)
                    tmp_gname = self.ambimat
                else:
                    n = glasscatalog.get_n(ele.data['material'][s], self.wl)
                    tmp_gname = ele.data['material'][s]
                if self.testmode:
                    if n == tmp_gname:
                        self.logfile.write(' ' + tmp_gname)
                self.data['type'].append(ele.data['type'][s])
                self.data['radius'].append(ele.data['radius'][s])
                self.data['asph'].append(ele.data['asph'][s])
                self.data['radius2'].append(ele.data['radius2'][s])
                self.data['asph2'].append(ele.data['asph2'][s])
                # self.data['freeform_coeff'].append(ele.data['freeform_coeff'][s])
                self.data['CT_sum'].append(CT_sum)
                self.data['lrad'].append(max_rCA)
                self.data['rCA'].append(ele.data['rCA'][s])
                self.data['rCA_short'].append(ele.data['rCA_short'][s])
                self.data['surf_rotation'].append(ele.data['surf_rotation'][s])
                self.data['CA_shape'].append('circle') # TODO
                self.data['surf_decenter'].append(ele.data['surf_decenter'][s])
                self.data['surf_tilt'].append(ele.data['surf_tilt'][s])
                self.data['n'].append(n)
                self.data['ismirror'].append(ele.ismirror or ele.data['coating'][s][0] == 'MIRROR')
                # coating
                coating_key = self._add_coating_data(ele.data['coating'][s])
                coating_entry = [coating_key, ele.data['coating'][s]]
                self.data['coating'].append(coating_entry)
                surfacecount = surfacecount + 1
                maxrad = max(maxrad, ele.data['rCA'][s])

        if self.imagemat:
            n = glasscatalog.get_n(self.imagemat, self.wl)
            self.data['n'][-1] = n
        self.isbuilt = True
        self.num_surfaces = surfacecount
        self.maximum_radius = maxrad

        # define detector quad
        self.set_detector()
        
        # add the sensor as a surface as well for handling sensor ghosts
        self.data['type'].append('flat')
        self.data['radius'].append(0)
        self.data['asph'].append([None, None])
        self.data['radius2'].append(0)
        self.data['asph2'].append([None, None])
        # self.data['freeform_coeff'].append([None])
        # self.data['CT'].append(None)
        self.data['CT_sum'].append(self.data['CT_sum'][-1] + self.detector['distance'])
        self.data['lrad'].append(max(self.detector['sizex'], self.detector['sizey']))
        self.data['rCA'].append(self.detector['sizex'])
        self.data['rCA_short'].append(self.detector['sizey'])
        self.data['surf_rotation'].append(0)
        self.data['CA_shape'].append('rect')
        self.data['surf_decenter'].append(None)
        self.data['surf_tilt'].append(None)
        self.data['n'].append(n)
        self.data['ismirror'].append(False)
        coating_key = self._add_coating_data(['FIXVALUE', 0.3])
        coating_entry = [coating_key, ['FIXVALUE', None, None]]
        self.data['coating'].append(coating_entry)

        # convert aperture positions to index
        # TODO: The other way round also?
        for i, ap in self.apertures.items():
            if ap['idx_surface'] is not None:
                # ignore apertures already defined by an index
                continue
            z_ap = ap['z_ap']
            if z_ap <= 0:
                self.apertures[i]['idx_surface'] = 0
                continue
            for j in range(1, self.num_surfaces):
                z_lens = self.data['CT_sum'][j]
                if z_ap > z_lens:
                    self.apertures[i]['idx_surface'] = j
                    

    #def add_aperture(self, aprad, zap, nblades):
    #    self.apertures.append([aprad, zap, nblades])

    #def define_sensor(self, zdet, detsizex, detsizey, relative=True):
    #    dz = 0
    #    if relative:
    #        dz = self.surfaces[-1][2]
    #    self.sensor = [detsizex, detsizey,zdet+dz]
    
    def set_detector(self):
        z = self.detector['distance'] + self.data['CT_sum'][-1]
        sx = self.detector['sizex']
        sy = self.detector['sizey']
        v1 = [-sx/2, -sy/2, z]
        v2 = [sx/2, -sy/2, z]
        v3 = [sx/2, sy/2, z]
        v4 = [-sx/2, sy/2, z]
        # pp = self.detector['pixelpitch']
        # v1 = [-sx*pp/2, -sy*pp/2, z]
        # v2 = [sx*pp/2, -sy*pp/2, z]
        # v3 = [sx*pp/2, sy*pp/2, z]
        # v4 = [-sx*pp/2, sy*pp/2, z]
        self.detector['Quad'] = [v1, v2, v3, v4]
        

    def change_wavelength(self, wl=0.546074):
        pass

    def reverse_system(self, build=True):
        pass

    def trace_surfaces_paraxial(self, y=1, u=0, axis=0):
        if not self.isbuilt:
            return None, None
        direction = 1
        for i in range(1, self.num_surfaces-1): # minus last surface and first dummy surface
            # surface
            if axis == 0:
                radius_use = 'radius'
            else:
                radius_use = 'radius2'
            r = self.data[radius_use][i]
            n1 = self.data['n'][i]
            n0 = self.data['n'][i-1]
            CT = self.data['CT_sum'][i+1] - self.data['CT_sum'][i]
            if self.data['ismirror'][i]:
                if r == 0:
                    pass
                else:
                    y, u = paraxial.reflect(y, u, r*direction)
                direction = direction*-1
            else:
                if r == 0:
                    phi = 0
                    #phi = -2./r
                else:
                    phi = (n1 - n0) / r
                y, u = paraxial.refract(y, u, n0, n1, phi)
            # thickness
            y, u = paraxial.transfer(y, u, CT*direction)
        # trace the final surface
        r = self.data[radius_use][i+1]
        n1 = self.data['n'][-1]
        n0 = self.data['n'][i]
        if self.data['ismirror'][i+1]:
            if r == 0:
                pass
            else:
                y, u = paraxial.reflect(y, u, r*direction)
            direction = direction*-1
        else:
            if r == 0:
                phi = 0
            else:
                phi = (n1 - n0) / r
            y, u = paraxial.refract(y, u, n0, n1, phi) # n_list[-1] ensures the ambient medium is always taken, so that an overly long list can be taken and the execution covered by n_elements
        
        return y, u
    
    def EFL_paraxial(self):
        if not self.isbuilt:
            return None
        y, u = self.trace_surfaces_paraxial(1, 0)
        EFL = paraxial.calc_EFL(1, u)
        return EFL

    def BFL_paraxial(self):
        if not self.isbuilt:
            return None
        y, u = self.trace_surfaces_paraxial(1, 0)
        BFL = paraxial.calc_BFL(y, u)
        return BFL
