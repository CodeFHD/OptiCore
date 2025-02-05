import numpy as np

import bpy
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from . import OBJECT_OT_add_lens, OBJECT_OT_add_sqlens#, OBJECT_OT_add_mirror, OBJECT_OT_add_CCretro

class OBJECT_OT_lens_system(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lenssystem"
    bl_label = "Lens System"
    bl_options = {'REGISTER', 'UNDO'}

    def create_system(system_name="Lens_System"):
        #create basic settings
        lens_dict = {}

        #add empty and name
        bpy.ops.object.add()
        bpy.context.active_object.name = system_name
        #if there was a duplicate swap around
        if not bpy.context.active_object.name == system_name:
            system_name = bpy.context.active_object.name
        #create custom property to say this is a lens system
        bpy.context.active_object["IS_OPTICORE_SYSTEM"] = "True"
        bpy.context.active_object["OPTICORE_SYSTEM_NAME"] = system_name

        #create individual lenses
        #and add name to empty custom property list



    def update_system():
        pass
        

    def execute(self, context):
        create_system(self, context)
        return {'FINISHED'}

def create_system(self, context):
    pass