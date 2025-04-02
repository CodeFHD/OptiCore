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

# from mathutils import Vector

import bpy
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from ..bl_materials import add_diffusematerial_cycles, add_diffusematerial_luxcore

def add_sensor(self, context, lx, ly, zsensor, thicksensor=False, sensorthickness=None):
    verts = []
    edges = []
    faces = []

    if thicksensor and sensorthickness is None:
        sensorthickness = max(lx, ly)/20

    # front face
    verts.append([0, -lx, -ly])
    verts.append([0, lx, -ly])
    verts.append([0, lx, ly])
    verts.append([0, -lx, ly])
    faces.append([0, 1, 2, 3])

    if thicksensor:
        # rear face
        verts.append([-sensorthickness, -lx, -ly])
        verts.append([-sensorthickness, lx, -ly])
        verts.append([-sensorthickness, lx, ly])
        verts.append([-sensorthickness, -lx, ly])
        faces.append([7, 6, 5, 4])
        # sides
        faces.append([0, 4, 5, 1]) # bottom
        faces.append([1, 5, 6, 2]) # right
        faces.append([2, 6, 7, 3]) # top
        faces.append([0, 3, 7, 4]) # left

    #create mesh from verts and faces
    mesh = bpy.data.meshes.new(name = 'OC_Sensor')
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)

    using_cycles = context.scene.render.engine == 'CYCLES'
    using_luxcore = context.scene.render.engine == 'LUXCORE'

    if using_cycles:
        materialname_sensor = add_diffusematerial_cycles(objectname='Sensor')
    elif using_luxcore:
        materialname_sensor = add_diffusematerial_luxcore(objectname='Sensor')
    if using_cycles or using_luxcore:
        material_sensor = bpy.data.materials[materialname_sensor]
        ob = bpy.context.active_object
        ob.data.materials.append(material_sensor)