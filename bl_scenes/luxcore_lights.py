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
import bpy

import numpy as np

from ..utils.geometry import get_rotmat_x, get_rotmat_z

def add_laser(location=(0,0,0), rotation=(0,0,0),
              lasersize=0.1):
    lightname = 'OC_Laser'
    light_data = bpy.data.lights.new(name=lightname, type='AREA')
    light_object = bpy.data.objects.new(name=lightname, object_data=light_data)
    bpy.context.collection.objects.link(light_object)

    light_object.data.luxcore.is_laser = True
    light_object.data.size = lasersize

    light_object.location = location

    # rotation is handled in YZX mode
    # Y should be 90 degree to align along the negative x-axis
    # Then Z describes the FOV angle
    # and X finally the azimuth
    light_object.rotation_mode = 'YZX'
    light_object.rotation_euler = rotation

    return light_object


def add_laser_array(HFOV=0, VFOV=0, numH=4, numV=4,
                    distance=20,
                    rp=[0,0,0], lasersize=0.1):
    created_lights = []
    if HFOV == 0:
        HFOV_list = [0]
    else:
        HFOV_list = np.linspace(-HFOV, HFOV, 2*numH+1, endpoint=True)*np.pi/180
    if VFOV == 0:
        VFOV_list = [0]
    else:
        VFOV_list = np.linspace(-VFOV, VFOV, 2*numV+1, endpoint=True)*np.pi/180

    for H in HFOV_list:
        for V in VFOV_list:
            # the following calculation is based on the assumption of
            # h = f*tan(theta), but
            # since f will cancel it is ommited here (or equivalently, assumed f=1)
            # step 1: determine rotation
            hx = np.tan(H)
            hy = np.tan(V)
            theta = np.arctan(np.sqrt(hx*hx + hy*hy))
            phi = np.arctan2(hy, hx)
            rotation = [phi, np.pi/2, theta]
            # step 2, rotate location around rp (rotation point)
            rp = np.array(rp)
            location = np.array([distance, 0, 0]) - rp # base location
            RZ = get_rotmat_z(theta)
            RX = get_rotmat_x(phi)
            location = np.matmul(RZ, location)
            location = np.matmul(RX, location)
            location = location + rp
            obj = add_laser(location, rotation, lasersize=lasersize)
            created_lights.append(obj.name)

    return created_lights