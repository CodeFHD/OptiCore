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

bl_info = {
    "name": "OptiCore",
    "author": "Johannes Hinrichs (CodeFHD)",
    "version": (2, 0),
    "blender": (4, 2, 5),
    "location": "View3D > Add > Mesh",
    "description": "Adds a new optical element",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh",
    }

import bpy
import numpy as np
# from bpy.types import Operator
# from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add # NOTE: have to import object_data_add here even though it is not used in this file, due to the way it is imported in other files. TODO: attempt to change it to absolute imports in those modules to clean this up.
from .elements import OBJECT_OT_add_lens, reset_lens_defaults, OBJECT_OT_add_sqlens, OBJECT_OT_add_mirror, OBJECT_OT_add_CCretro, OBJECT_OT_add_siemens
from .optomech import OBJECT_OT_add_table, OBJECT_OT_add_post
from .systems import OBJECT_OT_load_zmx, reset_loadzmx_faults#, OBJECT_OT_test_zmx

class OBJECT_OT_reset_performance_variables(bpy.types.Operator, AddObjectHelper):
    """This function resets some variables that may have caused performance issues"""
    bl_idname = "mesh.reset_performance_variables"
    bl_label = "Reset Performance Critical Settings"
    bl_options = {'REGISTER', 'UNDO'}
     
    def draw(self, context):
        pass

    def execute(self, context):
        reset_lens_defaults()
        reset_loadzmx_faults()
        return {'FINISHED'}

class OBJECT_MT_opticsmenu(bpy.types.Menu):
    bl_idname = 'mesh.opticsmenu'
    bl_label = 'OptiCore Optics'

    def draw(self, context):
        self.layout.operator(OBJECT_OT_add_lens.bl_idname)
        # self.layout.operator(OBJECT_OT_add_sqlens.bl_idname)
        self.layout.operator(OBJECT_OT_add_mirror.bl_idname)
        self.layout.operator(OBJECT_OT_add_CCretro.bl_idname)
        self.layout.operator(OBJECT_OT_add_siemens.bl_idname)
        self.layout.operator(OBJECT_OT_load_zmx.bl_idname)
        # self.layout.operator(OBJECT_OT_test_zmx.bl_idname)
        self.layout.operator(OBJECT_OT_reset_performance_variables.bl_idname)

class OBJECT_MT_optomechmenu(bpy.types.Menu):
    bl_idname = 'mesh.optomechmenu'
    bl_label = 'OptiCore Optomechanics'

    def draw(self, context):
        self.layout.operator(OBJECT_OT_add_table.bl_idname)
        self.layout.operator(OBJECT_OT_add_post.bl_idname)
    
def menu_func(self, context):
    self.layout.separator()
    self.layout.menu(OBJECT_MT_opticsmenu.bl_idname)
    self.layout.menu(OBJECT_MT_optomechmenu.bl_idname)

# Registration

def add_lens_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_lens.bl_idname,
        text="Add Lens",
        icon='PLUGIN')
"""
def add_lsqens_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_sqlens.bl_idname,
        text="Add Square Lens",
        icon='PLUGIN')
"""
def add_mirror_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_mirror.bl_idname,
        text="Add Mirror",
        icon='PLUGIN')
def add_table_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_table.bl_idname,
        text="Add Table",
        icon='PLUGIN')
def load_zmx_button(self, context):
    self.layout.operator(
        OBJECT_OT_load_zmx.bl_idname,
        text="Load lenses from .zmx file",
        icon='PLUGIN')
"""
def test_zmx_button(self, context):
    self.layout.operator(
        OBJECT_OT_test_zmx.bl_idname,
        text="Test zmx library",
        icon='PLUGIN')
"""
def reset_performance_button(self, context):
    self.layout.operator(
        OBJECT_OT_reset_performance_variables.bl_idname,
        text="Reset performance critical variables",
        icon='PLUGIN')

"""
# This allows you to right click on a button and link to the manual
def add_lens_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_lens", "editors/3dview/object"),
        )
    return url_manual_prefix, url_manual_mapping
def add_mirror_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_mirror", "editors/3dview/object"),
        )
    return url_manual_prefix, url_manual_mapping
"""

classes = (OBJECT_OT_reset_performance_variables,
           OBJECT_OT_add_lens, OBJECT_OT_add_mirror, OBJECT_OT_add_CCretro, # , OBJECT_OT_add_sqlens
           OBJECT_OT_add_siemens, OBJECT_OT_load_zmx, # OBJECT_OT_test_zmx,
          OBJECT_OT_add_table,OBJECT_OT_add_post)
menus = (OBJECT_MT_opticsmenu, OBJECT_MT_optomechmenu)

def register():
    for cla in classes+menus:
        print('REGISTERING',  cla)
        bpy.utils.register_class(cla)
    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)
    
def unregister():
    for cla in classes+menus:
        bpy.utils.unregister_class(cla)
    bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)

if __name__ == "__main__":
    register()
