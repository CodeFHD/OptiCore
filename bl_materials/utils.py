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

def clear_all_materials():
    # bpy.ops.outliner.orphans_purge() # TODO: This would be useful but might be messing with users scenes because it deletes too much. --> Limit to OC_ meshes and materials
    materials = bpy.data.materials
    for material in materials[:]: # use a copy so that removal doesn't mess with the iteration
        if material.users == 0 and material.name.startswith('OC_'):
            materials.remove(material)

def get_OC_material_name(material_name, engine):
    OC_material_name = 'OC_' + material_name + '_' + engine
    return OC_material_name


def check_OC_material(material_name, engine):
    OC_material_name = get_OC_material_name(material_name, engine)
    material_exists = OC_material_name in bpy.data.materials.keys()
    return material_exists, OC_material_name