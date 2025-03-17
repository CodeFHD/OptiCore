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

from .utils import *
from ..raytrace import glasscatalog


def add_glass_luxcore(glassname, n, luxcoreType='LuxCoreNodeMatGlass', update=False, reflection_color=[0, 0, 0]):
    OC_material_name = get_OC_material_name(glassname, 'luxcore')
    # Check if material already exists
    if OC_material_name in bpy.data.materials.keys() and not update:
        return bpy.data.materials.get(OC_material_name)
    elif update:
        material = bpy.data.materials.get(OC_material_name)
    else:
        material = bpy.data.materials.new(name=OC_material_name)
        material.use_nodes = True
    # set up a clear node tree
    luxcore_tree = bpy.data.node_groups.new(name='Nodes_' + OC_material_name, type="luxcore_material_nodes")
    material.luxcore.node_tree = luxcore_tree
    nodes = luxcore_tree.nodes
    nodes.clear()
    # define basic node components
    glass_node = nodes.new(type=luxcoreType)
    output_node = nodes.new(type='LuxCoreNodeMatOutput')
    # GUI layout
    glass_node.location = (0, 0)
    output_node.location = (300, 0)
    # set parameters
    glass_node.inputs["Transmission Color"].default_value = [1, 1, 1]
    glass_node.inputs["Reflection Color"].default_value = reflection_color
    glass_node.inputs["IOR"].default_value = n
    # connections
    luxcore_tree.links.new(glass_node.outputs["Material"], output_node.inputs["Material"])
    return material

def add_blackoutmaterial_luxcore(objectname='LensEdge'):
    OC_material_name = get_OC_material_name(objectname, 'luxcore')
    if OC_material_name in bpy.data.materials.keys():
        return OC_material_name# bpy.data.materials.get(OC_material_name)
    else:
        material = bpy.data.materials.new(name=OC_material_name)
        material.use_nodes = True
    # set up a clear node tree
    luxcore_tree = bpy.data.node_groups.new(name='Nodes_' + OC_material_name, type="luxcore_material_nodes")
    material.luxcore.node_tree = luxcore_tree
    nodes = luxcore_tree.nodes
    nodes.clear()
    # define basic node components
    matte_node = nodes.new(type='LuxCoreNodeMatMatte')
    output_node = nodes.new(type='LuxCoreNodeMatOutput')
    # GUI layout
    matte_node.location = (0, 0)
    output_node.location = (300, 0)
    # set parameters
    matte_node.inputs["Diffuse Color"].default_value = [0, 0, 0]
    material.diffuse_color = [0, 0, 0, 1]
    # connections
    luxcore_tree.links.new(matte_node.outputs["Material"], output_node.inputs["Material"])
    return OC_material_name

def add_diffusematerial_luxcore(objectname='LensDface', color=[1, 1, 1], viewportcolor=[1, 1, 1, 1]):
    OC_material_name = get_OC_material_name(objectname, 'luxcore')
    if OC_material_name in bpy.data.materials.keys():
        return OC_material_name# bpy.data.materials.get(OC_material_name)
    else:
        material = bpy.data.materials.new(name=OC_material_name)
        material.use_nodes = True
    # set up a clear node tree
    luxcore_tree = bpy.data.node_groups.new(name='Nodes_' + OC_material_name, type="luxcore_material_nodes")
    material.luxcore.node_tree = luxcore_tree
    nodes = luxcore_tree.nodes
    nodes.clear()
    # define basic node components
    diffuse_node = nodes.new(type='LuxCoreNodeMatMatte')
    output_node = nodes.new(type='LuxCoreNodeMatOutput')
    # GUI layout
    diffuse_node.location = (0, 0)
    output_node.location = (300, 0)
    # set parameters
    diffuse_node.inputs["Diffuse Color"].default_value = color
    material.diffuse_color = viewportcolor
    # connections
    luxcore_tree.links.new(diffuse_node.outputs["Material"], output_node.inputs["Material"])
    return OC_material_name

def glass_from_Element_luxcore(ele, wl, mat_refract_only=False):
    materials_bulk = []
    materials_interface = []
    num_glasses = len(ele.data['type']) - 1
    n_list = []

    luxcoreType = 'LuxCoreNodeMatGlass'
    if mat_refract_only:
        reflection_color = [0, 0, 0]
    else:
        reflection_color = [1, 1, 1]
    
    # create materials
    for i in range(num_glasses):
        # bulk material
        glassname = ele.data['material'][i]
        glassname_is_dummy = glassname.startswith('FIXVALUE_') or glassname.startswith('___BLANK')
        n = glasscatalog.get_n(glassname, wl)
        n_list.append(n)
        if glassname_is_dummy:
            glassname = f'Glass-{n:.5f}'
        if mat_refract_only:
            glassname = glassname + '_refraction'
        material_exists, OC_material_name = check_OC_material(glassname, 'luxcore')
        if not material_exists:
            _ = add_glass_luxcore(glassname, n, luxcoreType=luxcoreType, reflection_color=reflection_color)
        materials_bulk.append(OC_material_name)

    # interface materials, if present
    if num_glasses > 1:
        for i in range(num_glasses - 1):
            n_ratio = n_list[i+1]/n_list[i]
            glassname1 = ele.data['material'][i]
            glassname2 = ele.data['material'][i+1]
            glassname1_is_dummy = glassname1.startswith('FIXVALUE_') or glassname1.startswith('___BLANK')
            glassname2_is_dummy = glassname2.startswith('FIXVALUE_') or glassname2.startswith('___BLANK')
            if (glassname1_is_dummy or glassname2_is_dummy):
                glassname = f'GlassIF-{n_ratio:.5f}'
            else:
                glassname = f'GlassIF-{glassname1}-{glassname2}'
            material_exists, OC_material_name = check_OC_material(glassname, 'luxcore')
            if not material_exists:
                _ = add_glass_luxcore(glassname, n_ratio, luxcoreType=luxcoreType, reflection_color=reflection_color)
            materials_interface.append(OC_material_name)

    return materials_bulk, materials_interface