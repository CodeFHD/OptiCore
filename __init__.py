bl_info = {
    "name": "OptiCore",
    "author": "Johannes Hinrichs (CodeFHD)",
    "version": (1, 0),
    "blender": (2, 83, 0),
    "location": "View3D > Add > Mesh > Add Lens",
    "description": "Adds a new optical element",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh",
    }

import bpy
import numpy as np
from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from .elements import OBJECT_OT_add_lens, OBJECT_OT_add_sqlens, OBJECT_OT_add_mirror, OBJECT_OT_add_CCretro
from .optomech import OBJECT_OT_add_table, OBJECT_OT_add_post

class OBJECT_MT_opticsmenu(bpy.types.Menu):
    bl_idname = 'mesh.opticsmenu'
    bl_label = 'OptiCore Optics'

    def draw(self, context):
        self.layout.operator(OBJECT_OT_add_lens.bl_idname)
        self.layout.operator(OBJECT_OT_add_sqlens.bl_idname)
        self.layout.operator(OBJECT_OT_add_mirror.bl_idname)
        self.layout.operator(OBJECT_OT_add_CCretro.bl_idname)

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
def add_lsqens_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_sqlens.bl_idname,
        text="Add Square Lens",
        icon='PLUGIN')
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

classes = (OBJECT_OT_add_lens, OBJECT_OT_add_sqlens, OBJECT_OT_add_mirror, OBJECT_OT_add_CCretro,
          OBJECT_OT_add_table,OBJECT_OT_add_post)
menus = (OBJECT_MT_opticsmenu, OBJECT_MT_optomechmenu)

def register():
    for cla in classes+menus:
        bpy.utils.register_class(cla)
    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)
    
def unregister():
    for cla in classes+menus:
        bpy.utils.unregister_class(cla)
    bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)


#register, unregister = bpy.utils.register_classes_factory(classes)

if __name__ == "__main__":
    register()
