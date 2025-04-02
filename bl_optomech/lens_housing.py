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

import bpy
from bpy_extras.object_utils import object_data_add

from ..bl_materials import add_blackoutmaterial_cycles, add_blackoutmaterial_luxcore

def add_lenshousing_simple(self, context, lens, verts_outline, dz_outline,
                           CA_factor=0.99, housing_factor=1.1, thicksensor=False, sensorthickness=None,
                           dshape=False, outlinetype='max'):
    """
    This function creates a basic lens housing intended for quick raytracing setups.
    The basic logic is to create and aperture at every surface, scaled by CAfactor,
    and connecting them on the outside.
    Likewise, the sensor housing is modelled based on the sensor size and connected to the lens tube.
    """
    if dshape:
        # Not currently implemented for cross section lenses
        return
    if outlinetype not in ['max', 'tight']:
        print(f"[OC] Warning: Invalid value for outlinetype in add_lenshousing_simple(): '{outlinetype}'. Defaulting to 'max'")
        outlinetype = 'max'
    
    verts = []
    edges = []
    faces = []

    # create vertices at lenses and connect apertures
    offset = 0
    offset_list = [0]
    N_list = []
    i_surf = 0
    for vo, dz in zip(verts_outline, dz_outline):
        v, squarelens = vo
        N_this = len(v)
        N_list.append(N_this)
        if i_surf%2 == 0:
            # TODO: this is not ideal and not a generally good solution
            # However, a better solution might need to re-evaluate each surface twice. And this approach may acutally be OK.
            # Reason: the tube outside the lenses is most important, slight leakage effect at the baffles could be acceptable
            # If anything, solidyfying the tube might be more important regarding leaks
            dz_epsilon = 0.1
        else:
            dz_epsilon = -0.1
        v2 = np.array(v) # create a copy
        if outlinetype == 'max':
            if squarelens:
                thisrad = abs(v2[0, 1])
            else:
                thisrad = np.sqrt(np.sum(v2[0, 1:]*v2[0, 1:]))
            maxrad = np.nanmax(lens.data['lrad'][1:])
            outline_factor = maxrad/thisrad
        else:
            outline_factor = 1
        v2[:,1:] = v2[:,1:]*CA_factor
        v2[:,0] = v2[:,0] + dz + dz_epsilon
        v2 = [[v2[i,0], v2[i,1], v2[i,2]] for i in range(N_this)]
        verts = verts + v2
        v2 = np.array(v)
        v2[:,1:] = v2[:,1:]*housing_factor*outline_factor
        v2[:,0] = v2[:,0] + dz + dz_epsilon
        v2 = [[v2[i,0], v2[i,1], v2[i,2]] for i in range(N_this)]
        verts = verts + v2
        for j in range(N_this):
            fi1 = offset + j
            fi2 = offset + (j+1)%N_this
            fi3 = offset + (j+1)%N_this + N_this
            fi4 = offset + j + N_this 
            faces.append([fi4, fi3, fi2, fi1])
        offset += 2*N_this
        offset_list.append(offset)
        i_surf += 1

    # save some numbers for interfacing with sensor box
    v2 = np.array(v2)
    rad_last = np.sqrt(np.max(np.sum(v2[:,1:]**2, axis=1)))
    nVerts_tube = len(verts)
    z_tube_last = v2[0,0] + dz_epsilon # add dz_epsilon to avoid issues with bridge_edge_loops

    # connect lens tube
    for i in range(1, len(verts_outline)):
        if not N_list[i] == N_list[i-1]:
            # unequal sides connected later using bridge_edge_loops
            continue
        N_this = N_list[i]
        offset0 = offset_list[i-1] + N_this # index offset to the previous aperture incl. inner ring
        offset1 = offset_list[i] + N_this # index offset to this aperture incl. inner ring
        for j in range(N_this):
            fi1 = offset0 + j
            fi2 = offset0 + (j+1)%N_this
            fi3 = offset1 + (j+1)%N_this
            fi4 = offset1 + j
            faces.append([fi4, fi3, fi2, fi1])

    # create box around sensor
    length_lensbarrel = max(dz_outline) - min(dz_outline)
    # get sensor position and size
    # x/y extent
    len_sensor = max(lens.detector['sizex'], lens.detector['sizex'])
    len_sh = max(len_sensor/2, rad_last)*1.1
    # z-position. Assumption: Sersor is perpendicular to z-axis
    det_quad = np.array(lens.detector['Quad'])
    det_z = det_quad[0][2]
    if thicksensor:
        z_sh_back = -det_z - 2*sensorthickness
    else:
        z_sh_back = -det_z - 0.05*length_lensbarrel
    verts.append([z_tube_last, -len_sh, -len_sh])
    verts.append([z_tube_last, len_sh, -len_sh])
    verts.append([z_tube_last, len_sh, len_sh])
    verts.append([z_tube_last, -len_sh, len_sh])
    verts.append([z_sh_back, -len_sh, -len_sh])
    verts.append([z_sh_back, len_sh, -len_sh])
    verts.append([z_sh_back, len_sh, len_sh])
    verts.append([z_sh_back, -len_sh, len_sh])
    o = nVerts_tube # just for keeping the following lines short
    faces.append([o, o+1, o+5, o+4][::-1]) # sides
    faces.append([o+1, o+2, o+6, o+5][::-1])
    faces.append([o+2, o+3, o+7, o+6][::-1])
    faces.append([o+3, o, o+4, o+7][::-1])
    faces.append([o+4, o+5, o+6, o+7][::-1]) # rear wall
    
    # create object
    mesh = bpy.data.meshes.new(name = 'OC_LensHousing')
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    # bridge unconnected
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.context.tool_settings.mesh_select_mode = [True, False, False] # Vertex select mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    for i in range(1, len(verts_outline)):
        if N_list[i] == N_list[i-1]:
            # opposite case as in step "connect lens tube"
            continue
        N_before = N_list[i-1]
        N_this = N_list[i]
        offset0 = offset_list[i-1] + N_before # index offset to the previous aperture incl. inner ring
        offset1 = offset_list[i] + N_this # index offset to this aperture incl. inner ring
        for j in range(N_before):
            fi1 = offset0 + j
            mesh.vertices[fi1].select=True
        for j in range(N_this):
            fi3 = offset1 + (j+1)%N_this
            mesh.vertices[fi3].select=True

        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.bridge_edge_loops()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

    mesh = obj.data # get the updated mesh

    # create interface between sensor box and lens tube
    i = len(verts_outline) - 1 # should be preserved from above but bettzer be explicit
    N_this = N_list[i]
    offset1 = offset_list[i] + N_this # index offset to this aperture incl. inner ring
    for j in range(N_this):
        fi3 = offset1 + (j+1)%N_this
        mesh.vertices[fi3].select=True
    o = nVerts_tube
    for j in range(4):
        ei1 = o + j
        mesh.vertices[ei1].select=True
    
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    bpy.ops.mesh.bridge_edge_loops()
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

    # assign material
    using_cycles = context.scene.render.engine == 'CYCLES'
    using_luxcore = context.scene.render.engine == 'LUXCORE'
    
    if using_cycles:
        materialname_aperture = add_blackoutmaterial_cycles(objectname='Housing')
    elif using_luxcore:
        materialname_aperture = add_blackoutmaterial_luxcore(objectname='Housing')
    if using_cycles or using_luxcore:
        material_apertue = bpy.data.materials[materialname_aperture]
        ob = bpy.context.active_object
        ob.data.materials.append(material_apertue)