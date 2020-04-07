bl_info = {
    "name": "OptiCore",
    "author": "Johannes Hinrichs (CodeFHD)",
    "version": (2, 0),
    "blender": (2, 81, 0),
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
from . import elements as ele
from .elements.sqlens import OBJECT_OT_add_sqlens

class OBJECT_OT_add_lens(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "OptiCore Lens"
    bl_options = {'REGISTER', 'UNDO'}
    
    ltype1 : EnumProperty(
           name="Surface 1 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 1",
           #options={'HIDDEN'},
           )
    ltype2 : EnumProperty(
           name="Surface 2 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 2",
           #options={'HIDDEN'},
           )
    rad1 : FloatProperty(
           name="Surface 1 Radius",
           default = 12.,
           description="Radius of Curvature of Surface 1",
           unit = "LENGTH",
           )
    rad2 : FloatProperty(
           name="Surface 2 Radius",
           default = 24.,
           description="Radius of Curvature of Surface 2",
           unit = "LENGTH",
           )
    num1 : IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           min=3,
           )
    num2 : IntProperty(
           name="N2",
           default = 64,
           description="Number of angular vertices",
           min=3,
           )
    lensradius : FloatProperty(
           name="Lens Radius",
           default = 3.,
           description="Lens outer radius",
           unit = "LENGTH",
           )
    centerthickness : FloatProperty(
           name="Center Thickness",
           default = 1.,
           description="Center thickness of lens",
           unit = "LENGTH",
           )
    k : FloatProperty(
           name="k",
           default = 0.,
           description="Aspheric conical constant",
           )
    A : FloatVectorProperty(
           name="A",
           default = (0.,0.,0.),
           description="Aspheric correction coefficients",
           )
    k2 : FloatProperty(
           name="k2",
           default = 0.,
           description="Aspheric conical constant",
           )
    A2 : FloatVectorProperty(
           name="A2",
           default = (0.,0.,0.),
           description="Aspheric correction coefficients",
           )
    material_name : StringProperty(
            name="Material",
            default="",
           )
    shade_smooth : BoolProperty(
            name="Smooth Shading",
            default=True,
           )
    smooth_type : BoolProperty(
            name="Use Autosmooth (LuxCore v2.3)",
            default=True,
           )
    dshape : BoolProperty(
            name="D-Shaped Lens",
            default=False,
           )

    def draw(self, context):
        layout = self.layout
        # Location
        col = layout.column(align=True)
        col.label(text="Location")
        col.prop(self, 'location', text="")
        # Rotation
        col = layout.column(align=True)
        col.label(text="Rotation")
        col.prop(self, 'rotation', text="")
        scene = context.scene
        layout.prop(self, 'ltype1')
        layout.prop(self, 'rad1')
        layout.prop(self, 'ltype2')
        layout.prop(self, 'rad2')
        layout.prop(self, 'lensradius')
        layout.prop(self, 'centerthickness')
        layout.prop(self, 'num1')
        layout.prop(self, 'num2')
        if self.ltype1=='aspheric':
            layout.prop(self, 'k')
            layout.prop(self, 'A')
        if self.ltype2=='aspheric':
            layout.prop(self, 'k2')
            layout.prop(self, 'A2')
        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')
        layout.prop(self, 'dshape')

    def execute(self, context):
        ele.add_lens(self, context)
        return {'FINISHED'}

class OBJECT_OT_add_mirror(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_mirror"
    bl_label = "OptiCore Mirror"
    bl_options = {'REGISTER', 'UNDO'}
    
    mtype : EnumProperty(
           name="Surface Shape",
           items = {("parabolic","Parabolic",""),
                    ("spherical","Spherical","")},
           default = "parabolic",
           description="Shape of Mirror Surface",
           )
    opos : EnumProperty(
           name="Origin position",
           items = {("FP","Focal Point",""),
                    ("MC","Mirror Center","")},
           default = "FP",
           description="Position of the Mesh Origin w.r.t. optical properties",
           )
    rad : FloatProperty(
           name="Surface Radius",
           default = 12.,
           description="Radius of Curvature of Mirror Surface",
           unit = "LENGTH",
           )
    num1 : IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           min=3,
           )
    num2 : IntProperty(
           name="N2",
           default = 64,
           description="Nubmer of angular vertices",
           min=3,
           )
    mirrorradius : FloatProperty(
           name="Mirror Radius",
           default = 3.,
           description="Mirror outer radius",
           unit = "LENGTH",
           )
    centerthickness : FloatProperty(
           name="Back Thickness",
           default = 1.,
           description="Thickness at thinnest point",
           unit = "LENGTH",
           )
    theta : FloatProperty(
           name="Offset Angle",
           default = 0.,
           description="Offset angle for off-axis mirror",
           unit = "ROTATION",
           )
    material_name : StringProperty(
            name="Material",
            default="",
           )
    shade_smooth : BoolProperty(
            name="Smooth Shading",
            default=True,
           )
    smooth_type : BoolProperty(
            name="Use Autosmooth (LuxCore v2.3)",
            default=True,
           )
    cent_hole : BoolProperty(
            name="Central Hole",
            default=False,
           )
    hole_rad : FloatProperty(
           name="HoleRadius",
           default = 0.1,
           description="Radius of Curvature of Mirror Surface",
           min = 0.01,
           unit = "LENGTH",
           )

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        layout.prop(self, 'mtype')
        if self.mtype == 'parabolic':
            layout.prop(self, 'opos')
        layout.prop(self, 'rad')
        layout.prop(self, 'mirrorradius')
        layout.prop(self, 'centerthickness')
        if self.mtype == 'parabolic':
            layout.prop(self, 'theta')
        layout.prop(self, 'num1')
        layout.prop(self, 'num2')
        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        #layout.prop(self, 'material_name')
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')
        layout.prop(self, 'cent_hole')
        if self.cent_hole:
            layout.prop(self, 'hole_rad')

    def execute(self, context):

        ele.add_mirror(self, context)

        return {'FINISHED'}
    
def menu_func(self, context):
    self.layout.operator(OBJECT_OT_add_lens.bl_idname)
    self.layout.operator(OBJECT_OT_add_sqlens.bl_idname)
    self.layout.operator(OBJECT_OT_add_mirror.bl_idname)

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


classes = (OBJECT_OT_add_lens, OBJECT_OT_add_sqlens, OBJECT_OT_add_mirror)

def register():
    for cla in classes:
        bpy.utils.register_class(cla)
    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)
    
def unregister():
    for cla in classes:
        bpy.utils.unregister_class(cla)
    bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)


#register, unregister = bpy.utils.register_classes_factory(classes)

if __name__ == "__main__":
    register()
