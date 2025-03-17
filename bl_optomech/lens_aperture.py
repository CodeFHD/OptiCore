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

from mathutils import Vector

import bpy
from bpy_extras.object_utils import object_data_add

from ..bl_materials import add_blackoutmaterial_cycles, add_diffusematerial_cycles, add_blackoutmaterial_luxcore, add_diffusematerial_luxcore

def add_circular_aperture(self, context, radius_inner, radius_outer, N, dshape=False):
    verts = []
    edges = []
    faces = []
    
    maxb = 2*np.pi
    if dshape:
        maxb = np.pi*N/(N-1)
    
    # front face
    # create inner circle vertices
    for i in range(N):
        b = maxb*i/N
        verts.append(Vector((0, radius_inner*np.sin(b), radius_inner*np.cos(b))))
    # create outer circle vertices
    for i in range(N):
        b = maxb*i/N
        verts.append(Vector((0, radius_outer*np.sin(b), radius_outer*np.cos(b))))
    # back face
    # create inner circle vertices
    for i in range(N):
        b = maxb*i/N
        verts.append(Vector((-0.1, radius_inner*np.sin(b), radius_inner*np.cos(b))))
    # create outer circle vertices
    for i in range(N):
        b = maxb*i/N
        verts.append(Vector((-0.1, radius_outer*np.sin(b), radius_outer*np.cos(b))))
    
    # connect faces
    # front faces    
    for i in range(N - dshape):
        fi1 = i
        fi2 = (i+1)%N
        fi3 = i + N
        fi4 = (i+1)%N + N
        faces.append([fi1, fi2, fi4, fi3])
    # back faces    
    for i in range(N - dshape):
        fi1 = i + 2*N
        fi2 = (i+1)%N + 2*N
        fi3 = i + N + 2*N
        fi4 = (i+1)%N + N + 2*N
        faces.append([fi3, fi4, fi2, fi1])
    # side faces inner    
    for i in range(N - dshape):
        fi1 = i
        fi2 = (i+1)%N
        fi3 = i + 2*N
        fi4 = (i+1)%N + 2*N
        faces.append([fi3, fi4, fi2, fi1])
    # side faces outer    
    for i in range(N - dshape):
        fi1 = i + N
        fi2 = (i+1)%N + N
        fi3 = i + N + 2*N
        fi4 = (i+1)%N + N + 2*N
        faces.append([fi1, fi2, fi4, fi3])
    # d-shape
    if dshape:
        # top
        fi1 = 0
        fi2 = N
        fi3 = 2*N
        fi4 = 3*N
        faces.append([fi1, fi2, fi4, fi3])
        # bottom
        fi1 = N - 1
        fi2 = 2*N - 1
        fi3 = 3*N - 1
        fi4 = 4*N - 1
        faces.append([fi3, fi4, fi2, fi1])
    
    #create mesh from verts and faces
    mesh = bpy.data.meshes.new(name="Aperture")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)

    using_cycles = context.scene.render.engine == 'CYCLES'
    using_luxcore = context.scene.render.engine == 'LUXCORE'
    
    if using_cycles:
        materialname_aperture = add_blackoutmaterial_cycles(objectname='Aperture')
        if dshape: materialname_dface = add_diffusematerial_cycles(objectname='ApertureXSection')
    elif using_luxcore:
        materialname_aperture = add_blackoutmaterial_luxcore(objectname='Aperture')
        if dshape: materialname_dface = add_diffusematerial_luxcore(objectname='ApertureXSection')
    if using_cycles or using_luxcore:
        material_apertue = bpy.data.materials[materialname_aperture]
        ob = bpy.context.active_object
        ob.data.materials.append(material_apertue)
        if dshape:
            material_dface = bpy.data.materials[materialname_dface]
            ob.data.materials.append(material_dface)
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='DESELECT')
            sel_mode = bpy.context.tool_settings.mesh_select_mode
            bpy.context.tool_settings.mesh_select_mode = [False, False, True] # face select mode
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            mesh.polygons[-1].select = True
            mesh.polygons[-2].select = True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            obj.active_material_index = 1
            bpy.ops.object.material_slot_assign()
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            bpy.context.tool_settings.mesh_select_mode = [True, False, False] # reactivate vertex select mode