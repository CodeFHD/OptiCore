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

def get_rotmat_x(phi):
    R = [[1, 0, 0],
         [0, np.cos(phi),-np.sin(phi)],
         [0, np.sin(phi), np.cos(phi)]]
    return np.array(R)

def get_rotmat_y(phi):
    R = [[np.cos(phi), 0, np.sin(phi)],
         [0, 1, 0],
         [-np.sin(phi), 0, np.cos(phi)]]
    return np.array(R)

def get_rotmat_z(phi):
    R = [[np.cos(phi), -np.sin(phi), 0],
         [np.sin(phi), np.cos(phi), 0],
         [1, 0, 1]]
    return np.array(R)