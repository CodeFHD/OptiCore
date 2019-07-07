bl_info = {
    "name": "OptiCore",
    "author": "Johannes Hinrichs (CodeHD)",
    "version": (1, 0),
    "blender": (2, 79, 0),
    "location": "View3D > Add > Mesh > Add Lens",
    "description": "Adds a new optical element",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh",
    }

import bpy
import numpy as np
from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from . import elements as ele

class OBJECT_OT_add_lens(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}
    
    rad1 = FloatProperty(
           name="Surface 1 Radius",
           default = 12.,
           description="Radius of Curvature of Surface 1",
           )
    rad2 = FloatProperty(
           name="Surface 2 Radius",
           default = 24.,
           description="Radius of Curvature of Surface 2",
           )
    num1 = IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           )
    num2 = IntProperty(
           name="N2",
           default = 64,
           description="Number of angular vertices",
           )
    lensradius = FloatProperty(
           name="Lens Radius",
           default = 3.,
           description="Lens outer radius",
           )
    centerthickness = FloatProperty(
           name="Center Thickness",
           default = 1.,
           description="Center thickness of lens",
           )
    material_name = StringProperty(
            name="Material",
            default="",
           )
    shade_smooth = BoolProperty(
            name="Smooth Shading",
            default=True,
           )
    split_edge = BoolProperty(
            name="Edge Split",
            default=True,
           )

    def execute(self, context):

        ele.add_lens(self, context)

        return {'FINISHED'}

class OBJECT_OT_add_mirror(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_mirror"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}
    
    mtype = EnumProperty(
           name="Surface Shape",
           items = {("parabolic","Parabolic",""),
                    ("spherical","Spherical","")},
           default = "parabolic",
           description="Radius of Curvature of Mirror Surface",
           )
    opos = EnumProperty(
           name="Origin position",
           items = {("FP","Focal Point",""),
                    ("MC","Mirror Center","")},
           default = "FP",
           description="Position of the Mesh Origin w.r.t. optical properties",
           )
    rad = FloatProperty(
           name="Surface 1 Radius",
           default = 12.,
           description="Radius of Curvature of Mirror Surface",
           )
    num1 = IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           )
    num2 = IntProperty(
           name="N2",
           default = 64,
           description="Nubmer of angular vertices",
           )
    mirrorradius = FloatProperty(
           name="Lens Radius",
           default = 3.,
           description="Mirror outer radius",
           )
    centerthickness = FloatProperty(
           name="Thickness",
           default = 1.,
           description="Thickness at thinnest point",
           )
    theta = FloatProperty(
           name="Offset Angle",
           default = 0.,
           description="Offset angle for off-axis mirror",
           )

    def execute(self, context):

        ele.add_mirror(self, context)

        return {'FINISHED'}

# Registration

def add_lens_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_lens.bl_idname,
        text="Add Lens",
        icon='PLUGIN')
def add_mirror_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_mirror.bl_idname,
        text="Add Mirror",
        icon='PLUGIN')


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


def register():
    bpy.utils.register_class(OBJECT_OT_add_mirror)
    bpy.utils.register_manual_map(add_mirror_manual_map)
    bpy.types.INFO_MT_mesh_add.append(add_mirror_button)
    bpy.utils.register_class(OBJECT_OT_add_lens)
    bpy.utils.register_manual_map(add_lens_manual_map)
    bpy.types.INFO_MT_mesh_add.append(add_lens_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_mirror)
    bpy.utils.unregister_manual_map(add_mirror_manual_map)
    bpy.types.INFO_MT_mesh_add.remove(add_mirror_button)
    bpy.utils.unregister_class(OBJECT_OT_add_lens)
    bpy.utils.unregister_manual_map(add_lens_manual_map)
    bpy.types.INFO_MT_mesh_add.remove(add_lens_button)


if __name__ == "__main__":
    register()
