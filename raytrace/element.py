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

"""
This file implements an optical Element class.
An Element is equivalent to what is usally called a group in a lens system, i.e. a singlet, a cemented doublet, etc.
Another example is a single mirror in a reflective telescope.

The surfaces in an Element are defined relative to the optical axis of the first surface.
In other words, the Element is defined in its local coordinates and does not store information about global 3D-positon and orientation.

The separate class Lenssystem forms a complete system based on a list of Elements and their geometric arrangement. 
"""

ALLOWED_SURFTYPES = ['flat',                                                 # flat shapes
    'spherical',    'conical',          'polynominal',      'aspheric',      # axial symmetry
    'cylindrical',  'conicylindrical',  'polycylindrical',  'acylindrical',  # cylinder symmtry
    'toric',        'conitoric',        'polytoric',        'atoric']        # toric shape

class Element():
    def __init__(self, name='NONAME'):
        self.name = name
        self.data = {}
        
        """lists with length n_surfaces"""
        ## surface geometry - optical
        self.data['type'] = [] # Surface type for specific functions
        self.data['lenstype'] = [] # For the OptiCore Lens class, stores its shape class, e.g. cylindricalX
        self.data['radius'] = [] # Spherical surface radius [r0, r1, ...]. 
        self.data['asph'] = [] # Apsheric coefficients [[k, a4, a6,...], ...].
        self.data['radius2'] = [] # secondary radius list for toric lenses
        self.data['asph2'] = [] # secondary aspheric list for toric lenses
        # self.data['freeform_coeff'] = [] # coefficients for freeform-equation

        ## surface geometry - mechanical
        self.data['lrad'] = [] # Mechanical Lens radius. Radius for round lenses, side half-length for square lens. List of two half-lengths for rectangular lenses
        self.data['rCA'] = [] # Radius of optically ground surface. Same logic as for lrad. Must be smaller or equal to lrad. Difference will be created as flat annulus
        self.data['rCA_short'] = [] # Short axis of rCA, used in case of e.g. cylinder lens
        self.data['surf_rotation'] = [] # angle by which surface is rotated around optical axis. E.g. for cylinder lenses.

        ## imperfections
        self.data['surf_decenter'] = [] # Surface lateral position error [x, y] in local coordinates
        self.data['surf_tilt'] = [] # Surface tilt error [X, Y] in local coordinates

        ## other surface properties
        self.data['coating'] = [] # coatings [type, datasource, dataindex]
        # self.data['roughness_model'] = [] # for scattering
        # self.data['roughness_coeffs']
        
        # """lists with length n_surfaces - 1"""
        self.data['CT'] = [] # Center thickness
        self.data['material'] = [] # Glass type specifiers
        
        # TEMPLATE
        # self.data[''] = []

        # global element parameters
        self.outline_shape = 'circular' # alternative options: square, rectangular
        
        # status flags
        self.ismirror = False # This is specifically in case of a single, discrete mirror with one surface
        self.isdummy = False # dummy lenses for placing e.g. an aperture
        self.isreverse = False # indicates if the element has been reversed from its original state
        self.direction = 1 # indicates to the Blender object generation in which direction the element is facing
        # self.orientation = [0, 0, 0] # angles in spherical coordiantes [latitude, longitude, axis roll] to specify how this Element shall be oriented in 3D-space. Order of application: 1. Roll aorund z-axis, 2. latitude around x-axis, 3. longitude around y-axis.

    def add_surface(self, surftype, lenstype='rotational', radius=None, asph=[None, None], radius2=None, asph2=[None, None],
                    coating=[None, None, None], rCA=None, rCA_short = None, lrad=None, CT=None, material=None,
                    surf_decenter=None, surf_tilt=None, surf_rotation=0):
        """
        function to append appropriately to the data structure
        """
        if not surftype in ALLOWED_SURFTYPES:
            print('ERROR: Element.add_surface received unknown surftype:', surftype)
            return

        if len(self.data['radius']) != len(self.data['CT']):
            print('ERROR in Element.add_surface: Mismatch of data structure length. Possibly trying to append to completed element?')
            return

        if rCA_short is None: rCA_short = rCA

        # parameters independent of surftype
        self.data['type'].append(surftype)
        self.data['lenstype'].append(lenstype)
        self.data['lrad'].append(lrad)
        self.data['rCA'].append(rCA)
        self.data['rCA_short'].append(rCA_short)
        self.data['surf_decenter'].append(surf_decenter)
        self.data['surf_tilt'].append(surf_tilt)
        self.data['coating'].append(coating)
        if CT is not None:
            self.data['CT'].append(CT)
            self.data['material'].append(material)

        if surftype == 'flat':
            self.data['radius'].append(0)
            self.data['asph'].append([None, None])
            self.data['radius2'].append(None)
            self.data['asph2'].append([None, None])
            self.data['surf_rotation'].append(0)
        elif surftype in ['spherical', 'conical', 'polynominal', 'aspheric']:
            self.data['radius'].append(radius)
            self.data['asph'].append(asph)
            self.data['radius2'].append(None)
            self.data['asph2'].append([None, None])
            self.data['surf_rotation'].append(0)
        elif surftype in ['cylindrical', 'conicylindrical', 'polycylindrical', 'acylindrical',
                          'toric', 'conitoric', 'polytoric', 'atoric']:
            """
            Also for cylindric lenses, second variables have to be specified explicitly because cylinder might be along x- or y-axis.
            TODO: replace with surf_rotation parameter.
            """
            self.data['radius'].append(radius)
            self.data['asph'].append(asph)
            self.data['radius2'].append(radius2)
            self.data['asph2'].append(asph2)
            self.data['surf_rotation'].append(surf_rotation)

    def reverse(self):
        # TODO: This needs to be overhauled to account for latest set of parameters
        for key in self.data:
            if key in ['radius', 'radius2']:
                self.data[key] = [-1.*r for r in self.data[key][::-1]]
            elif key in ['asph', 'asph2'] and self.data[key] != []:
                self.data[key] = [[e[0]] + [-1.*e[i] for i in range(1,len(e))] for e in self.data[key][::-1]] # conical constant remains unchanged
            elif key == 'surf_tilt' and self.data[key] != []:
                self.data[key] = [[-1.*e[0], -1.*e[1]] for e in self.data[key][::-1]]
            #elif key == 'glass':
            #    self.data[key] = [self.data[key][:-1:-1]] + [self.data[key][-1]]
            else:
                self.data[key] = self.data[key][::-1]
        self.isreverse =  not self.isreverse
        
    def decenter(self):
        pass