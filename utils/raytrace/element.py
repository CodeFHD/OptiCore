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

class Element():
    def __init__(self, name='NONAME'):
        self.name = name
        self.data = {}
        
        # lists with length n_surfaces
        self.data['type'] = [] # Surface type for specific functions
        self.data['radius'] = [] # Spherical surface radius
        self.data['asph'] = [] # Apsheric coefficients [[k, a4, a6, ...]]
        self.data['coating'] = [] # coatings
        # self.data['lrad'] = [] # Mechanical Lens radius
        self.data['rCA'] = [] # Clear aperture radius
        self.data['surf_decenter'] = [] # Surface lateral position error [[x, y]] in local coordinates
        self.data['surf_tilt'] = [] # Surface tilt error [[X, Y]] in local coordinates
        
        # lists with length n_surfaces - 1
        self.data['CT'] = [] # Center thickness
        self.data['material'] = [] # Glass type specifiers
        
        # status flags
        self.ismirror = False # This is specifically in case of a single, discrete mirror with one surface
        self.isdummy = False # dummy lenses for placing e.g. an aperture
        self.isreverse = False # indicates if the element has been reversed from its original state
        self.direction = 1 # indicates to the Blender object generation in which direction the element is facing
        # self.orientation = [0, 0, 0] # angles in spherical coordiantes [latitude, longitude, axis roll] to specify how this Element shall be oriented in 3D-space. Order of application: 1. Roll aorund z-axis, 2. latitude around x-axis, 3. longitude around y-axis.

    def reverse(self):
        for key in self.data:
            if key == 'radius':
                self.data[key] = [-1.*r for r in self.data[key][::-1]]
            elif key == 'asph' and self.data[key] != []:
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
    