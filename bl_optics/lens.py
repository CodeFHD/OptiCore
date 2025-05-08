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
from bpy.types import SEQUENCER_PT_sequencer_overlay_strips
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils
from ..raytrace import paraxial
from ..raytrace.element import Element
from ..raytrace.lenssystem import Lenssystem
from ..raytrace.trace_sequential import exec_trace
from ..bl_raytrace.trace_scene import trace_to_scene
from ..raytrace import rayfan
from ..utils.check_surface import surftype_Lens2Element
from ..surface.surfaceutils import get_N1_sqsurface

def reset_lens_defaults():
    bpy.context.window_manager.operator_properties_last("mesh.add_lens").num1 = 32
    bpy.context.window_manager.operator_properties_last("mesh.add_lens").num2 = 64
    bpy.context.window_manager.operator_properties_last("mesh.add_lens").nrays = 11
    bpy.context.window_manager.operator_properties_last("mesh.add_lens").addrayfan = False
    bpy.context.window_manager.operator_properties_last("mesh.add_lens").tracetoscene = False

class OBJECT_OT_add_lens(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "Lens"
    bl_options = {'REGISTER', 'UNDO'}
    
    display_option : EnumProperty(
        name="Settings",
        items = (("default","Most Important",""),
                 ("optical","Optical Design",""),
                 ("geometry","Model Geometry",""),
                 ("rayfan","Ray Tracing",""),
                 ("paraxial","Paraxial Optics",""),),
        default = "default",
        description="Select which settings are displayed in this window.",
           )
    """Optical layout"""
    makedoublet : EnumProperty(
        name="Lens configuration",
        items = (("1", "Singlet", ""),
                 ("2", "Doublet", ""),
                 ("3", "Triplet", "")),
        default="1",
        )
    squarelens : BoolProperty(
            name="Square Lens",
            default=False,
           )
    ltype_options = (("flat", "Flat", ""),
                    ("rotational","Standard",""),
                    ("cylindrical","Cylindrical",""),
                    #("cylindricY","Cylindrical (Y)",""),
                    ("toric", "Toric", ""))
    ltype1 : EnumProperty(
           name="",
           items = ltype_options,
           default = "rotational",
           description="Shape of Surface 1",
           #options={'HIDDEN'},
           )
    ltype2 : EnumProperty(
           name="",
           items = ltype_options,
           default = "rotational",
           description="Shape of Surface 2",
           #options={'HIDDEN'},
           )
    ltype3 : EnumProperty(
           name="",
           items = ltype_options,
           default = "rotational",
           description="Shape of Surface 3",
           #options={'HIDDEN'},
           )
    ltype4 : EnumProperty(
           name="",
           items = ltype_options,
           default = "rotational",
           description="Shape of Surface 4",
           #options={'HIDDEN'},
           )
    """Surface parameters (optical)"""
    # primary radius
    RY1 : FloatProperty(
           name="",
           default = 61.47,
           description="Radius of Curvature",
           unit = "LENGTH",
           )
    RY2 : FloatProperty(
           name="",
           default = -44.64,
           description="Radius of Curvature of Surface 2",
           unit = "LENGTH",
           )
    RY3 : FloatProperty(
           name="",
           default = -129.94,
           description="Radius of Curvature of Surface 3",
           unit = "LENGTH",
           )
    RY4 : FloatProperty(
           name="",
           default = -126.,
           description="Radius of Curvature of Surface 4",
           unit = "LENGTH",
           )
    # secondary radius
    RZ1 : FloatProperty(
           name="",
           default = 61.47,
           description="Second radius of curvature (for cylindric/toric lenses)",
           unit = "LENGTH",
           )
    RZ2 : FloatProperty(
           name="",
           default = 61.47,
           description="Second radius of curvature (for cylindric/toric lenses)",
           unit = "LENGTH",
           )
    RZ3 : FloatProperty(
           name="",
           default = 61.47,
           description="Second radius of curvature (for cylindric/toric lenses)",
           unit = "LENGTH",
           )
    RZ4 : FloatProperty(
           name="",
           default = 61.47,
           description="Second radius of curvature (for cylindric/toric lenses)",
           unit = "LENGTH",
           )
    # surface rotation
    surfrot1 : FloatProperty(
           name="",
           default = 0,
           min = -2*np.pi,
           max = 2*np.pi,
           description="Rotation of surface around x-axis",
           unit = "ROTATION",
           )
    surfrot2 : FloatProperty(
           name="",
           default = 0,
           min = -2*np.pi,
           max = 2*np.pi,
           description="Rotation of surface around x-axis",
           unit = "ROTATION",
           )
    surfrot3 : FloatProperty(
           name="",
           default = 0,
           min = -2*np.pi,
           max = 2*np.pi,
           description="Rotation of surface around x-axis",
           unit = "ROTATION",
           )
    surfrot4 : FloatProperty(
           name="",
           default = 0,
           min = -2*np.pi,
           max = 2*np.pi,
           description="Rotation of surface around x-axis",
           unit = "ROTATION",
           )
    # aspheric coefficients
    #THIS DIDN'T WORK
    #ASPDEG : IntProperty(
    #       name='ASPDEG',
    #       default = 3,
    #       description = "Maximum Degree of Aspheric Coefficients to be Displayed",
    #       min = 3,
    #       max = 12,
    #       )
    ASPDEG = 3
    k1 : FloatProperty(
           name="",
           default = 0.,
           description="Conic constant",
           )
    A1 : FloatVectorProperty(
           name="",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    k2 : FloatProperty(
           name="",
           default = 0.,
           description="Aspheric conical constant",
           )
    A2 : FloatVectorProperty(
           name="",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    k3 : FloatProperty(
           name="",
           default = 0.,
           description="Aspheric conical constant",
           )
    A3 : FloatVectorProperty(
           name="",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    k4 : FloatProperty(
           name="",
           default = 0.,
           description="Aspheric conical constant",
           )
    A4 : FloatVectorProperty(
           name="",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    """Surface parameters (mechanical)"""
    lensradius : FloatProperty(
           name="Lens Radius",
           default = 12.5,
           description="Lens outer radius",
           unit = "LENGTH",
           )
    flangerad1 : FloatProperty(
           name="",
           default = 0.,
           description="Flange Width",
           unit = "LENGTH",
           min = 0,
           )
    flangerad2 : FloatProperty(
           name="",
           default = 0.,
           description="Flange Width Surface 2",
           unit = "LENGTH",
           min = 0,
           )
    flangerad3 : FloatProperty(
           name="",
           default = 0.,
           description="Flange Width Surface 3",
           unit = "LENGTH",
           min = 0,
           )
    flangerad4 : FloatProperty(
           name="",
           default = 0.,
           description="Flange Width Surface 4",
           unit = "LENGTH",
           min = 0,
           )
    centerthickness1 : FloatProperty(
           name="Center Thickness 1",
           default = 6.,
           description="Center thickness of lens",
           unit = "LENGTH",
           )
    centerthickness2 : FloatProperty(
           name="Center Thickness 2",
           default = 2.5,
           description="Center thickness of lens segment 2",
           unit = "LENGTH",
           )
    centerthickness3 : FloatProperty(
           name="Center Thickness 3",
           default = 4.,
           description="Center thickness of lens segment 3",
           unit = "LENGTH",
           )
    """System analysis (output only)"""
    def get_flen(self):
        return self.flen_intern
    flen_intern : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Paraxial EFL of the lens",
           )
    flen : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Paraxial EFL of the lens",
           get=get_flen,
           )
    def get_rmsspotsize(self):
        return self.rmsspotsize_intern
    rmsspotsize_intern : FloatProperty(
           name="RMS spot size",
           default = 0.,
           description="RMS spot size of the ray trace in the sensor plane",
           )
    rmsspotsize : FloatProperty(
           name="RMS spot size",
           default = 0.,
           description="RMS spot size of the ray trace in the sensor plane",
           get=get_rmsspotsize,
           )
    """Sequential ray tracing"""
    wavelength : StringProperty(
            name="Wavelengths",
            default="e;F';C'",
            description="Wavelengths for calculations and ray tracing. Values are separated by semicolon. The first value is taken as the primary wavelength. Values can either be entered in nanometers or using Fraunhofer line letters ",
           )
    ior1 : FloatProperty(
           name="IOR 1",
           default = 1.5168, # N-BK7 Schott
           description="Index of Refraction of first lens section. Not transferred to material, only for paraxial infromations.",
           )
    ior2 : FloatProperty(
           name="IOR 2",
           default = 1.67271, # N-SF5 Schott
           description="Index of Refraction of second lens section. Not transferred to material, only for paraxial infromations.",
           )
    ior3 : FloatProperty(
           name="IOR 3",
           default = 1.62005, # N-F2 Schott
           description="Index of Refraction of third lens section. Not transferred to material, only for paraxial infromations.",
           )
    addrayfan : BoolProperty(
            name="Add Ray Fan",
            default=False,
           )
    zdet : FloatProperty(
           name="Detector Distance",
           default = 0.0, # 48.0,
           description="Distance offset to the plane where rays are traced to. When this is set to 0, the paraxial BFL is used.",
           )
    nrays : IntProperty(
           name="Nrays",
           default = 11,
           description="Number of rays for ray fan",
           min=3,
           )
    fantype : EnumProperty(
           name="Ray Fan Type",
           items = (("2D","2D",""),
                    ("2D_finite","2D - finite",""),
                    ("3D_tri","3D tris",""),
                    ("3D_tri_finite","3D tris - finite",""),
                    ("3D_rings", "3D rings", ""),
                    ("3D_rings_finite", "3D rings - finite", ""),
                    ("3D_random","3D random",""),
                    ("3D_random_finite","3D random - finite",""),
                    ("3D_random_sun","3D random - sun",""),),
           default = "2D",
           description="Ray Fan Type",
           #options={'HIDDEN'},
           )
    fandist : FloatProperty(
           name="Ray Fan Origin",
           default = 20.,
           description="Distance where ray fan originates.",
           )
    fanangle1 : FloatProperty(
           name="Fan FoV Angle",
           default = 0.,
           description="Field-of-View angle of the ray fan.",
           min = -np.pi/2,
           max = np.pi/2,
           unit = "ROTATION",
           )
    fanangle2 : FloatProperty(
           name="Fan Azimuth Angle",
           default = np.pi/2,
           description="Azimuth angle of the ray fan.",
           min = 0,
           max = 2*np.pi,
           unit = "ROTATION",
           )
    fanangle3 : FloatProperty(
           name="Fan Roll Angle",
           default = np.pi/2.,
           description="Roll angle of the ray fan.",
           min = 0.,
           max = np.pi,
           unit = "ROTATION",
           )
    fandiam : FloatProperty(
           name="Ray Fan diameter factor",
           default = 0.9,
           description="Relative diameter of ray fan relaitve to first surface.",
           min = 0.01,
           max = 1.0,
           )
    tracetoscene : BoolProperty(
            name="Trace Ray Fan to Scene",
            default=False,
           )
    ghost_order : StringProperty(
            name="Ghost sequence",
            description="Specify two indices, seprated by a comma, between which surfaces a ghost reflection shall be traced. First surface has index 1.",
            default="",
           )
    """Optimization parameters"""
    autofocus : BoolProperty(
            name="RMS-Focus",
            description='Attempt to autofocus by ray tracing and evaluating the RMS spot size',
            default=False,
           )
    autosolve : BoolProperty(
            name="Automatically optimize one parameter",
            default=False,
           )
    solve_for : EnumProperty(
        name="Autosolver",
        items = (("r1","Surface 1 Radius",""),
                 ("r2","Surface 2 Radius",""),
                 ("r3","Surface 3 Radius",""),
                 ("r4","Surface 4 Radius",""),
                 ("n1","Glass 1 index",""),
                 ("n2","Glass 2 index",""),
                 ("n3","Glass 3 index",""),),
        default = "r1",
        description="Select which parameter to automatically optimize.",
           )
    flen_goal : FloatProperty(
           name="Focal Length target",
           default = 100.,
           description="Paraxial EFL of the lens used in autosolver.",
           min=1.,
           )
    """Mesh options"""
    material_name1 : StringProperty(
            name="Material1",
            description="Material Surface 1",
            default="",
           )
    material_name2 : StringProperty(
            name="Material2",
            description="Material Surface 2",
            default="",
           )
    material_name3 : StringProperty(
            name="Material3",
            description="Material Surface 3",
            default="",
           )
    material_name4 : StringProperty(
            name="Material4",
            description="Material Surface 4",
            default="",
           )
    num1 : IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           min=6,
           )
    num2 : IntProperty(
           name="N2",
           default = 64,
           description="Number of angular vertices",
           min=6,
           )
    shade_smooth : BoolProperty(
            name="Smooth Shading",
            default=True,
           )
    smooth_type : BoolProperty(
            name="Use Custom Normals",
            default=True,
           )
    dshape : BoolProperty(
            name="Cross-section Model",
            default=False,
           )
    display_edit : BoolProperty(
            name="Display Edit Mode",
            default=False,
           )

    def draw(self, context):
        md = self.makedoublet
        disp = self.display_option
        scene = context.scene
        row = self.layout
        col = row.column(align=True)
        if disp == 'default':
            col.prop(self, 'flen')
            col.prop(self, 'rmsspotsize')
            col.prop(self, 'lensradius')
            col.prop(self, 'RY1')
            col.prop(self, 'RY2')
            col.prop(self, 'RY3')
            col.prop(self, 'RY4')
            col.prop(self, 'centerthickness1')
            col.prop(self, 'centerthickness2')
            col.prop(self, 'centerthickness3')
            col.prop(self, 'ior1')
            col.prop(self, 'ior2')
            col.prop(self, 'ior3')
        if disp == 'optical':
            col.label(text="Lens Parameters")
            col.prop(self, 'flen')
            col.prop(self, 'rmsspotsize')
            col.prop(self, 'lensradius')
            #
            box = col.box()
            lrow = box.row()
            lrow.alignment="CENTER"
            lrow.label(text="Surface 1")
            brow = box.split(factor=0.2)
            bcol1 = brow.column(align=True)
            bcol1.label(text="Shape")
            bcol1.label(text="Rotation")
            bcol1.label(text="Radius")
            bcol1.label(text="\u03BA")
            bcol1.label(text="A4")
            bcol1.label(text="A6")
            bcol1.label(text="A8")
            bcol1.label(text="Flange")
            bcol2 = brow.column(align=True)
            bcol2.prop(self, 'ltype1')
            bcol2.prop(self, 'surfrot1', emboss=self.ltype1 in ['cylindrical', 'toric'])
            radrow = bcol2.row()
            radrow.prop(self, 'RY1')
            radrow.prop(self, 'RZ1', emboss=self.ltype1 in ['toric'])
            bcol2.prop(self, 'k1')
            bcol2.prop(self, 'A1')
            bcol2.prop(self, 'flangerad1')
            #
            col.prop(self, 'centerthickness1')
            #
            box = col.box()
            lrow = box.row()
            lrow.alignment="CENTER"
            lrow.label(text="Surface 2")
            brow = box.split(factor=0.2)
            bcol1 = brow.column(align=True)
            bcol1.label(text="Shape")
            bcol1.label(text="Rotation")
            bcol1.label(text="Radius")
            bcol1.label(text="\u03BA")
            bcol1.label(text="A4")
            bcol1.label(text="A6")
            bcol1.label(text="A8")
            bcol1.label(text="Flange")
            bcol2 = brow.column(align=True)
            bcol2.prop(self, 'ltype2')
            bcol2.prop(self, 'surfrot2', emboss=self.ltype2 in ['cylindrical', 'toric'])
            radrow = bcol2.row()
            radrow.prop(self, 'RY2')
            radrow.prop(self, 'RZ2', emboss=self.ltype2 in ['toric'])
            bcol2.prop(self, 'k2')
            bcol2.prop(self, 'A2')
            bcol2.prop(self, 'flangerad2')
            #
            col.prop(self, 'centerthickness2')
            #
            box = col.box()
            lrow = box.row()
            lrow.alignment="CENTER"
            lrow.label(text="Surface 3")
            brow = box.split(factor=0.2)
            bcol1 = brow.column(align=True)
            bcol1.label(text="Shape")
            bcol1.label(text="Rotation")
            bcol1.label(text="Radius")
            bcol1.label(text="\u03BA")
            bcol1.label(text="A4")
            bcol1.label(text="A6")
            bcol1.label(text="A8")
            bcol1.label(text="Flange")
            bcol2 = brow.column(align=True)
            bcol2.prop(self, 'ltype3')
            bcol2.prop(self, 'surfrot3', emboss=self.ltype3 in ['cylindrical', 'toric'])
            radrow = bcol2.row()
            radrow.prop(self, 'RY3')
            radrow.prop(self, 'RZ3', emboss=self.ltype3 in ['toric'])
            bcol2.prop(self, 'k3')
            bcol2.prop(self, 'A3')
            bcol2.prop(self, 'flangerad3')
            #
            col.prop(self, 'centerthickness3')
            #
            box = col.box()
            lrow = box.row()
            lrow.alignment="CENTER"
            lrow.label(text="Surface 4")
            brow = box.split(factor=0.2)
            bcol1 = brow.column(align=True)
            bcol1.label(text="Shape")
            bcol1.label(text="Rotation")
            bcol1.label(text="Radius")
            bcol1.label(text="\u03BA")
            bcol1.label(text="A4")
            bcol1.label(text="A6")
            bcol1.label(text="A8")
            bcol1.label(text="Flange")
            bcol2 = brow.column(align=True)
            bcol2.prop(self, 'ltype4')
            bcol2.prop(self, 'surfrot4', emboss=self.ltype4 in ['cylindrical', 'toric'])
            radrow = bcol2.row()
            radrow.prop(self, 'RY4')
            radrow.prop(self, 'RZ4', emboss=self.ltype4 in ['toric'])
            bcol2.prop(self, 'k4')
            bcol2.prop(self, 'A4')
            bcol2.prop(self, 'flangerad4')
        if disp == 'geometry':
            # Location
            col.label(text="Location")
            col.prop(self, 'location', text="")
            # Rotation
            col.label(text="Rotation")
            col.prop(self, 'rotation', text="")
            col.label(text="Modelling Parameters")
            col.prop(self, 'num1', emboss=not self.squarelens)
            col.prop(self, 'num2')
            col.prop_search(self, "material_name1", bpy.data, "materials", icon="NONE")
            col.prop_search(self, "material_name2", bpy.data, "materials", icon="NONE")
            col.prop_search(self, "material_name3", bpy.data, "materials", icon="NONE")
            col.prop_search(self, "material_name4", bpy.data, "materials", icon="NONE")
            col.prop(self, 'shade_smooth')
        if disp == 'rayfan':
            col.prop(self, 'flen')
            col.prop(self, 'rmsspotsize')
            col.prop(self, 'zdet')
            col.prop(self, 'nrays')
            col.prop(self, 'fantype')
            col.prop(self, 'fandist')
            col.prop(self, 'fanangle1')
            col.prop(self, 'fanangle2')
            col.prop(self, 'fanangle3')
            col.prop(self, 'fandiam')
            col.prop(self, 'tracetoscene')
            col.prop(self, 'autofocus')
            col.prop(self, 'ghost_order')
        if disp == 'paraxial':
            col.label(text="Optical Parameters")
            # col.prop(self, 'wavelength')
            col.prop(self, 'ior1')
            col.prop(self, 'ior2')
            col.prop(self, 'ior3')
            col.prop(self, 'flen')
            col.prop(self, 'rmsspotsize')
            
        # global config, always displayed
        col.prop(self, 'display_edit')
        col.prop(self, 'dshape')
        col.prop(self, 'squarelens')
        col.prop(self, 'addrayfan')
        col.prop(self, 'makedoublet')
        col.prop(self, 'display_option')
            
    def execute(self, context):
        add_lens(self, context)
        if self.addrayfan:
        #     utils.trace_rays(self, context)
            rmsspotsize = add_rayfan(self, context)
            self.rmsspotsize_intern = rmsspotsize
        return {'FINISHED'}
    
    @classmethod
    def reset_properties(cls):
        for prop, default_value in cls._defaults.items():
            setattr(cls, prop, default_value)
        # Reset each property to its default value
        # cls.centerthickness1 = cls.__annotations__['centerthickness1'].default
   



def get_default_paramdict_lens():
    """
    This functions returns a parameter dictionary filled with some default values.
    This can be used as a template for external calls to add_lens()
    when it is desired to only use a subset of the possible varaibles
    """
    paramdict = {}
    paramdict['squarelens'] = False
    paramdict['ltype1'] = "rotational"
    paramdict['ltype2'] = "rotational"
    paramdict['ltype3'] = "rotational"
    paramdict['ltype4'] = "rotational"
    paramdict['RY1'] = 12.
    paramdict['RY2'] = 24.
    paramdict['RY3'] = 24.
    paramdict['RY4'] = 24.
    paramdict['RZ1'] = 12.
    paramdict['RZ2'] = 24.
    paramdict['RZ3'] = 24.
    paramdict['RZ4'] = 24.
    paramdict['surfrot1'] = 0.
    paramdict['surfrot2'] = 0.
    paramdict['surfrot3'] = 0.
    paramdict['surfrot4'] = 0.
    paramdict['num1'] = 32
    paramdict['num2'] = 64
    paramdict['lensradius'] = 3.
    paramdict['flangerad1'] = 0.
    paramdict['flangerad2'] = 0.
    paramdict['flangerad3'] = 0.
    paramdict['flangerad4'] = 0.
    paramdict['centerthickness1'] = 1.
    paramdict['centerthickness2'] = 1.
    paramdict['centerthickness3'] = 1.
    paramdict['k1'] = 0.
    paramdict['A1'] = list([0. for i in range(3)])
    paramdict['k2'] = 0.
    paramdict['A2'] = list([0. for i in range(3)])
    paramdict['k3'] = 0.
    paramdict['A3'] = list([0. for i in range(3)])
    paramdict['k4'] = 0.
    paramdict['A4'] = list([0. for i in range(3)])
    paramdict['makedoublet'] = "1"
    paramdict['dshape'] = False
    paramdict['material_name1'] = ''
    paramdict['material_name2'] = ''
    paramdict['material_name3'] = ''
    paramdict['material_name4'] = ''
    paramdict['material_edge'] = ''
    paramdict['material_dface'] = ''
    paramdict['shade_smooth'] = True
    paramdict['smooth_type'] = True
    paramdict['ior1'] = 1.5
    paramdict['ior2'] = 1.6
    paramdict['ior3'] = 1.55
    paramdict['display_edit'] = False
    return paramdict

def add_lens(self, context, paramdict=None):
    # extract parameters
    if paramdict is None:
        squarelens = self.squarelens
        ltype1 = self.ltype1
        ltype2 = self.ltype2
        ltype3 = self.ltype3
        ltype4 = self.ltype4
        RY1 = self.RY1
        RY2 = self.RY2
        RY3 = self.RY3
        RY4 = self.RY4
        RZ1 = self.RZ1
        RZ2 = self.RZ2
        RZ3 = self.RZ3
        RZ4 = self.RZ4
        surfrot1 = self.surfrot1
        surfrot2 = self.surfrot2
        surfrot3 = self.surfrot3
        surfrot4 = self.surfrot4
        N1 = self.num1
        N2 = self.num2
        lrad = self.lensradius
        flrad1 = self.flangerad1
        flrad2 = self.flangerad2
        flrad3 = self.flangerad3
        flrad4 = self.flangerad4
        CT1 = self.centerthickness1
        CT2 = self.centerthickness2
        CT3 = self.centerthickness3
        k1 = self.k1
        A1 = self.A1
        k2 = self.k2
        A2 = self.A2
        k3 = self.k3
        A3 = self.A3
        k4 = self.k4
        A4 = self.A4
        md = self.makedoublet
        dshape = self.dshape
        material_name1 = self.material_name1
        material_name2 = self.material_name2
        material_name3 = self.material_name3
        material_name4 = self.material_name4
        material_edge = ''
        material_dface = ''
        shade_smooth = self.shade_smooth
        smooth_type = self.smooth_type
        ior1 =  self.ior1
        ior2 =  self.ior2
        ior3 =  self.ior3
        display_edit = self.display_edit
    else:
        squarelens = paramdict['squarelens']
        ltype1 = paramdict['ltype1']
        ltype2 = paramdict['ltype2']
        ltype3 = paramdict['ltype3']
        ltype4 = paramdict['ltype4']
        RY1 = paramdict['RY1']
        RY2 = paramdict['RY2']
        RY3 = paramdict['RY3']
        RY4 = paramdict['RY4']
        RZ1 = paramdict['RZ1']
        RZ2 = paramdict['RZ2']
        RZ3 = paramdict['RZ3']
        RZ4 = paramdict['RZ4']
        surfrot1 = paramdict['surfrot1']
        surfrot2 = paramdict['surfrot2']
        surfrot3 = paramdict['surfrot3']
        surfrot4 = paramdict['surfrot4']
        N1 = paramdict['num1']
        N2 = paramdict['num2']
        lrad = paramdict['lensradius']
        flrad1 = paramdict['flangerad1']
        flrad2 = paramdict['flangerad2']
        flrad3 = paramdict['flangerad3']
        flrad4 = paramdict['flangerad4']
        CT1 = paramdict['centerthickness1']
        CT2 = paramdict['centerthickness2']
        CT3 = paramdict['centerthickness3']
        k1 = paramdict['k1']
        A1 = paramdict['A1']
        k2 = paramdict['k2']
        A2 = paramdict['A2']
        k3 = paramdict['k3']
        A3 = paramdict['A3']
        k4 = paramdict['k4']
        A4 = paramdict['A4']
        md = paramdict['makedoublet']
        dshape = paramdict['dshape']
        material_name1 = paramdict['material_name1']
        material_name2 = paramdict['material_name2']
        material_name3 = paramdict['material_name3']
        material_name4 = paramdict['material_name4']
        material_edge = paramdict['material_edge']
        material_dface = paramdict['material_dface']
        shade_smooth = paramdict['shade_smooth']
        smooth_type = paramdict['smooth_type']
        ior1 =  paramdict['ior1']
        ior2 =  paramdict['ior2']
        ior3 =  paramdict['ior3']
        display_edit = paramdict['display_edit']

    """
    Verification of the geometry
    and resolving of degenerate cases
    """
    ret_dummy = utils.check_surface(lrad, flrad1, RY1, k1, A1, RZ1, None, None, ltype1, squarelens)
    lrad1, hasfl1, flrad1, ltype1, surf_subtype_X1, surf_subtype_Y1, RY1, k1, A1, dSurfrot1 = ret_dummy
    ret_dummy = utils.check_surface(lrad, flrad2, RY2, k2, A2, RZ2, None, None, ltype2, squarelens)
    lrad2, hasfl2, flrad2, ltype2, surf_subtype_X2, surf_subtype_Y2, RY2, k2, A2, dSurfrot2 = ret_dummy
    ret_dummy = utils.check_surface(lrad, flrad3, RY3, k3, A3, RZ3, None, None, ltype3, squarelens)
    lrad3, hasfl3, flrad3, ltype3, surf_subtype_X3, surf_subtype_Y3, RY3, k3, A3, dSurfrot3 = ret_dummy
    ret_dummy = utils.check_surface(lrad, flrad4, RY4, k4, A4, RZ4, None, None, ltype4, squarelens)
    lrad4, hasfl4, flrad4, ltype4, surf_subtype_X4, surf_subtype_Y4, RY4, k4, A4, dSurfrot4 = ret_dummy

    surfrot1 = surfrot1 + dSurfrot1
    surfrot2 = surfrot2 + dSurfrot2
    surfrot3 = surfrot3 + dSurfrot3
    surfrot4 = surfrot4 + dSurfrot4

    if material_edge == '':
        material_edge = material_name1
    if material_dface == '':
        material_dface = material_edge

    """
    prepare some variables
    """
    num_sides = int(md)       
    num_surfaces = int(md) + 1
    if squarelens:
        if N2%2 == 0:
            N2 = N2 + 1
        N1 = get_N1_sqsurface(N2, dshape)
    
    # collect into lists for looping
    ltype_list = [ltype1, ltype2, ltype3, ltype4]
    lsubtype_list = [surf_subtype_X1, surf_subtype_X2, surf_subtype_X3, surf_subtype_X4]
    srad_list = [RY1, RY2, RY3, RY4]
    RZ_list = [RZ1, RZ2, RZ3, RZ4]
    surfrot_list = [surfrot1, surfrot2, surfrot3, surfrot4]
    k_list = [k1, k2, k3, k4]
    A_list = [A1, A2, A3, A4]
    lrad_list = [lrad1, lrad2, lrad3, lrad4]
    hasfl = [hasfl1, hasfl2, hasfl3, hasfl4]
    CT_list = [0, CT1, CT2, CT3]
    
    # for keeping up with counts
    # named segment because sufacesd and sides alternate S1, S2, E1, S3, E2, ...
    # initialize with zero in order to facilitate taking differences starting from surface 0
    nVerts_at_segment = [0]
    nLoops_at_segment = [0]
    nFaces_at_segment = [0]

    # For lens housing
    verts_outline = []
    
    """
    Generating the vertices and faces
    """

    verts, faces, normals, N_inside_sq_prev, vo = sfc.add_surface(ltype1, surf_subtype_X1, squarelens, N1, N2, lrad1,
                                                          RY1, k=k1, A=A1,
                                                          RY2=RZ1, k2=None, A2=None,
                                                          surf_rotation=surfrot1,
                                                          dshape=dshape, lrad_ext=lrad)

    verts_outline.append([vo, squarelens])

    nVerts1 = len(verts)
    nFacs1 = len(faces)
    nVerts_at_segment.append(len(verts))
    nFaces_at_segment.append(len(faces))
    
    # all sides are equal so this is generic and gets added below
    if squarelens:
        pass
        # normalsside = sfc.get_sqringnormals(N1, N2, dshape=dshape)
    else:
        normalsside = sfc.get_ringnormals(N2, dshape=dshape)
        normalsside = [[t[2], t[0], t[1]] for t in normalsside]
    nFacSide = N2
    
    # add further surfaces
    CT_sum = 0
    srad_prev = RY1
    nVerts_tot = nVerts1
    nVerts_prev = nVerts1 # N2
    for i in range(1, int(md)+1):
        # get parameters for this surface
        # they are named with index 0 to avoid confusion with globals without number where that exists, e.g. lrad
        srad0 = srad_list[i]
        srad0Z = RZ_list[i]
        lrad0 = lrad_list[i]
        CT_sum = CT_sum + CT_list[i]
        ltype0 = ltype_list[i]
        surf_subtype0 = lsubtype_list[i]
        k0 = k_list[i]
        A0 = A_list[i]
        surfrot0 = surfrot_list[i]
        # if k0 == 0 and np.all(np.array(A0) == 0):
        #     ltype0 = 'spherical'
        
        # get the new surface and append

        dvert, dfac, dnormals, N_inside_sq_here, vo = sfc.add_surface(ltype0, surf_subtype0, squarelens, N1, N2, lrad0,
                                                          srad0, k=k0, A=A0,
                                                          RY2=srad0Z, k2=None, A2=None,
                                                          surf_rotation=surfrot0,
                                                          zadd=CT_sum, nVerts=nVerts_tot,
                                                          dshape=dshape, lrad_ext=lrad)

        # flip normals for last surface
        if i==int(md):
            dfac = [df[::-1] for df in dfac]
            dnormals = [(-n[0], -n[1], -n[2]) for n in dnormals]

        verts = verts + dvert
        faces = faces + dfac
        normals = normals + dnormals
        verts_outline.append([vo, squarelens])
        
        nVerts_this = len(dvert)
        nVerts_at_segment.append(len(verts))
        nFaces_at_segment.append(len(faces))
        
        # add the side
        if squarelens:
            # left - the first row of verts of each surface
            face_tmp1 = []
            face_tmp2 = []
            for j in range(N2):
                face_tmp1.append(nVerts_tot + j)
                face_tmp2.append(nVerts_tot + j - nVerts_prev)
            faces.append(face_tmp1[::-1] + face_tmp2)
            # right - the last row of verts of each surface
            face_tmp1 = []
            face_tmp2 = []
            for j in range(N2):
                face_tmp1.append(len(verts) - j - 1)
                face_tmp2.append(len(verts) - j - 1 - nVerts_this)
            faces.append(face_tmp1[::-1] + face_tmp2)
            # top
            face_tmp1 = []
            face_tmp2 = []
            for j in range(N1-1):
                face_tmp1.append(nVerts_tot + (N2-1) + N_inside_sq_here*j)
                face_tmp2.append(nVerts_tot + (N2-1) + N_inside_sq_prev*j - nVerts_prev)
            # last face needs special treatment, because row is always full even for flat surface
            face_tmp1.append(nVerts_tot + (N2-1) + N_inside_sq_here*j + N2)
            face_tmp2.append(nVerts_tot + (N2-1) + N_inside_sq_prev*j + N2 - nVerts_prev)
            faces.append(face_tmp1[::-1] + face_tmp2)
            # bottom
            face_tmp1 = []
            face_tmp2 = []
            # first face needs special treatment, because row is always full even for flat surface
            face_tmp1.append(nVerts_tot)
            face_tmp2.append(nVerts_tot - nVerts_prev)
            for j in range(N1-1):
                face_tmp1.append(nVerts_tot + N2 + N_inside_sq_here*j)
                face_tmp2.append(nVerts_tot + N2 + N_inside_sq_prev*j - nVerts_prev)
            faces.append(face_tmp1 + face_tmp2[::-1])

        else:
            for j in range(N2-dshape):
                fi1 = len(verts) - N2 + (j+1)%(N2)
                fi2 = len(verts) - N2 + j
                fi3 = fi2 - nVerts_this
                fi4 = fi1 - nVerts_this
                faces.append([fi4,fi3,fi2,fi1])
            if dshape:
                if srad_prev==0:
                    ptss1 = [nVerts_tot - 1, nVerts_tot - N2]
                else:
                    p0 = nVerts_tot - nVerts_prev 
                    ptss1 = [p0 + (x+1)*N2 for x in range(N1)][::-1] + [p0] + [p0 + 1 + x*N2 for x in range(N1)]
                if srad0==0:
                    ptss2 = [len(verts)-N2, len(verts)-1]
                else:
                    p0 = len(verts) - nVerts_this
                    ptss2 = [p0 + 1 + x*N2 for x in range(N1)][::-1] + [p0] + [p0 + (x+1)*N2 for x in range(N1)]
                dface = ptss1 + ptss2
                dface = dface[::-1]
                faces.append(dface)
        
        # update needed carry-over-parameters after this surface
        nVerts_prev = nVerts_this
        nVerts_tot = len(verts) # the (total) length of verts so far
        nVerts_at_segment.append(len(verts))
        nFaces_at_segment.append(len(faces))
        srad_prev = 1.*srad0
        if squarelens:
            N_inside_sq_prev = N_inside_sq_here
    
    #create mesh from verts and faces
    del dvert
    del dfac
    edges = [] # edges are not explicitly created here so we pass an empty list
    mesh = bpy.data.meshes.new(name = 'OC_Lens')
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    """
    Adding Blender Materials
    """
    #assign material(s)
    hasmat1 = material_name1 in bpy.data.materials
    hasmat2 = material_name2 in bpy.data.materials
    hasmat3 = material_name3 in bpy.data.materials
    hasmat4 = material_name4 in bpy.data.materials
    hasmat_edge = material_edge in bpy.data.materials
    hasmat_dface = material_dface in bpy.data.materials
    hasmat_any = hasmat1 or hasmat2 or hasmat3 or hasmat4 # not testing edge an dface because they are supplemental
    
    # initially change to edit mode
    if hasmat_any:
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [False, False, True] # face select mode
    
    # apend the materials to the objects material list
    if hasmat1:
        mat1 = bpy.data.materials[material_name1]
        if obj.data.materials.find(material_name1) < 0: # technically the check is unneccessary for material1, but added for consistency
            obj.data.materials.append(mat1) 
    if hasmat2:
        mat2 = bpy.data.materials[material_name2]
        if obj.data.materials.find(material_name2) < 0:
            obj.data.materials.append(mat2)
    if hasmat3:
        mat3 = bpy.data.materials[material_name3]
        if obj.data.materials.find(material_name3) < 0:
            obj.data.materials.append(mat3)
    if hasmat4:
        mat4 = bpy.data.materials[material_name4]
        if obj.data.materials.find(material_name4) < 0:
            obj.data.materials.append(mat4)
    if hasmat_edge:
        mat5 = bpy.data.materials[material_edge]
        if obj.data.materials.find(material_edge) < 0:
            obj.data.materials.append(mat5)
    if hasmat_dface and dshape:
        mat6 = bpy.data.materials[material_dface]
        if obj.data.materials.find(material_dface) < 0:
            obj.data.materials.append(mat6)
        
    # S1
    if hasmat1:    
        # deselect all
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        # select the relevant faces
        nfaces = nFaces_at_segment[1] - nFaces_at_segment[0]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False) # this needs to happen in object mode for some reason
        for i in range(nfaces):
            mesh.polygons[i].select=True
        # assign the material
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        obj.active_material_index = obj.data.materials.find(material_name1)
        bpy.ops.object.material_slot_assign()
    
    # S2
    if hasmat2:    
        # deselect all
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        # select the relevant faces
        nfaces = nFaces_at_segment[2] - nFaces_at_segment[1]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        for i in range(nfaces):
            mesh.polygons[i + nFaces_at_segment[1]].select=True
        # assign the material
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        obj.active_material_index = obj.data.materials.find(material_name2)
        bpy.ops.object.material_slot_assign()
    
    # S3
    if hasmat3 and md in ['2', '3']:  
        # deselect all
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        # select the relevant faces
        nfaces = nFaces_at_segment[4] - nFaces_at_segment[3]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        for i in range(nfaces):
            mesh.polygons[i + nFaces_at_segment[3]].select=True
        # assign the material
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        obj.active_material_index = obj.data.materials.find(material_name3)
        bpy.ops.object.material_slot_assign()
    
    # S4
    if hasmat4 and md in ['3']:  
        # deselect all
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        # select the relevant faces
        nfaces = nFaces_at_segment[6] - nFaces_at_segment[5]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        for i in range(nfaces):
            mesh.polygons[i + nFaces_at_segment[5]].select=True
        # assign the material
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        obj.active_material_index = obj.data.materials.find(material_name4)
        bpy.ops.object.material_slot_assign()
    
    # side
    if hasmat_edge:
        # deselect all
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        # select the relevant faces
        # side 1
        nfaces = nFaces_at_segment[3] - nFaces_at_segment[2] - dshape
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        for i in range(nfaces):
            mesh.polygons[i + nFaces_at_segment[2]].select=True
        # side 2
        if md in ['2', '3']:
            nfaces = nFaces_at_segment[5] - nFaces_at_segment[4] - dshape
            for i in range(nfaces):
                mesh.polygons[i + nFaces_at_segment[4]].select=True
        # side 3
        if md in ['3']:
            nfaces = nFaces_at_segment[5] - nFaces_at_segment[4] - dshape
            for i in range(nfaces):
                mesh.polygons[i + nFaces_at_segment[4]].select=True
            for i in range(nfaces):
                mesh.polygons[i + nFaces_at_segment[5]].select=True
        # assign the material
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        obj.active_material_index = obj.data.materials.find(material_edge)
        bpy.ops.object.material_slot_assign()
        # flanges
        nfaces = N2 - dshape
        if hasfl1 and not squarelens:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            for i in range(nfaces):
                mesh.polygons[nFaces_at_segment[1] - i - 1].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            obj.active_material_index = obj.data.materials.find(material_edge)
            bpy.ops.object.material_slot_assign()
        if hasfl2 and not squarelens:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            for i in range(nfaces):
                mesh.polygons[nFaces_at_segment[2] - i - 1].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            obj.active_material_index = obj.data.materials.find(material_edge)
            bpy.ops.object.material_slot_assign()
        if hasfl3 and not squarelens and md in ['2', '3']:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            for i in range(nfaces):
                mesh.polygons[nFaces_at_segment[4] - i - 1].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            obj.active_material_index = obj.data.materials.find(material_edge)
            bpy.ops.object.material_slot_assign()
        if hasfl4 and not squarelens and md in ['3']:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            for i in range(nfaces):
                mesh.polygons[nFaces_at_segment[6] - i - 1].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            obj.active_material_index = obj.data.materials.find(material_edge)
            bpy.ops.object.material_slot_assign()
    
    # d-shape
    if hasmat_dface and dshape:  
        # deselect all
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        # select the relevant faces
        # side 1
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        mesh.polygons[nFaces_at_segment[3] - 1].select=True
        # side 2
        if md in ['2', '3']:
            mesh.polygons[nFaces_at_segment[5] - 1].select=True
        # side 3
        if md in ['3']:
            mesh.polygons[nFaces_at_segment[7] - 1].select=True
        # assign the material
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        obj.active_material_index = obj.data.materials.find(material_dface)
        bpy.ops.object.material_slot_assign()
    
    """
    Creating split edge normals for smooth shading
    """

    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    """Case 1: Round lens"""
    if shade_smooth and not squarelens:
        bpy.ops.object.shade_smooth()
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
        bpy.ops.mesh.customdata_custom_splitnormals_add()
        cn_list = []
        i_segment = 1 # start at 1 because of 0 offset
        n_loops = 0

        for i_surf in range(num_surfaces):
            surf_subtype0 = lsubtype_list[i_surf]
            ltype0 = ltype_list[i_surf]
            # surface
            if surf_subtype0 == 'flat':
                nloops1 = N2 # when radius == 0, the face is flat circle, only one loop line going around
                nloops1_fl = N2 # flange has no impact for a flat face
            elif ltype0 in ['rotational']:
                nloops1 = (N2-dshape)*(3 + 4*(N1-1)) # center of lens has triangular faces, hence the 3. Rest quads, hence the 4. N-1 because of the innter triangle section.
                nloops1_fl = nloops1 - hasfl[i_surf]*4*(N2-dshape) # 2*hasfl[i_surf]*(3 + 4*(N1-1)+1)
            else:
                nloops1 = (N2-dshape)*(3 + 3*2*(N1-1)) # In this case, all triangular faces except flange
                nloops1_fl = nloops1 - hasfl[i_surf]*4*(N2-dshape) # 2*hasfl[i_surf]*(3 + 4*(N1-1)+1)
            for i in range(nloops1_fl):
                vi = mesh.loops[i + n_loops].vertex_index
                cn_list.append(normals[vi])
            if hasfl[i_surf]:
                for i in range(4*(N2-dshape)): # (2*(3 + 4*(N1-1))+2):
                    vi = mesh.loops[nloops1_fl + n_loops + i].vertex_index
                    cn_list.append(normals[nVerts_at_segment[i_segment]-1]) # all annulus vertices have the same pointing
            i_segment = i_segment + 1
            n_loops = n_loops + nloops1
            
            # no sides at the first face
            if i_surf == 0:
                continue
        
            # sides
            n_loop_side = 4*(N2 - dshape)
            nVerts = nVerts_at_segment[i_segment - 2]
            for i in range(n_loop_side):
                vi = mesh.loops[i + n_loops].vertex_index
                if vi < nVerts:
                    vi = (vi - nVerts + N2)%N2
                else:
                    vi = (vi - nVerts_at_segment[i_segment - 1])%N2
                cn_list.append(normalsside[vi])
            n_loops = n_loops + n_loop_side
        
            # D-shape
            if dshape:
                if srad_list[i_surf-1] == 0:
                    n_side1 = 1
                else:
                    n_side1 = 2*N1
                if srad_list[i_surf] == 0:
                    n_side2 = 1
                else:
                    n_side2 = 2*N1
                for i in range(n_side1 + n_side2 + 2):
                    cn_list.append((0,-1,0))
                n_loops = n_loops + n_side1 + n_side2 + 2
            
            i_segment = i_segment + 1

        mesh.normals_split_custom_set(cn_list)

    """Case 2: Square lens"""
    if shade_smooth and squarelens:
        bpy.ops.object.shade_smooth()
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
        bpy.ops.mesh.customdata_custom_splitnormals_add()
        cn_list = []
        i_segment = 1 # start at 1 because of 0 offset
        n_loops = 0

        for i_surf in range(num_surfaces):
            surf_subtype0 = lsubtype_list[i_surf]
            ltype0 = ltype_list[i_surf]
            # surface
            if surf_subtype0 == 'flat':
                nloops1 = 2*(N1 - 1) + 2*(N2 - 1) # when radius == 0, the face is flat square, edge has N1-1 or N2-1 loops
                nloops1_fl = nloops1 # flange has no impact for a flat face
            else:
                # all faces are tris for the sqlens
                nloops1 = 3*2*(N1 - 1)*(N2 - 1) # In this case, all triangular faces
                nloops1_fl = nloops1 # Needs to change when flanges implemented for sqlens
            for i in range(nloops1_fl):
                vi = mesh.loops[i + n_loops].vertex_index
                cn_list.append(normals[vi])
            i_segment = i_segment + 1
            n_loops = n_loops + nloops1

            # no sides at the first face
            if i_surf == 0:
                continue

            # side
            # explanation for n_loop_side: left/right (forst term) and top/bottom are identical so factor two each
            # then each side (factor 2 for a total of four) has N2-1 loops at the surface
            # plus the two lines connecting front and back face, cancelling out 2*(N2-1) + 2 = 2*N2
            # left face is identical with and withouot dshape
            n_loop_side =  + 4*N2 + 4*N1
            #Not working with normalsside here because it is easier to separate by side since they each have equal normal
            # left
            for i in range(2*N2):
                cn_list.append([0, -1, 0])
            # right
            for i in range(2*N2):
                cn_list.append([0, 1, 0])
            # top
            for i in range(2*N1):
                cn_list.append([0, 0, 1])
            # bottom
            for i in range(2*N1):
                cn_list.append([0, 0, -1])
            n_loops = n_loops + n_loop_side

            # D-shape
            # dshape does not need special handling for sqlens,
            # included in terms above because of necessity

            i_segment = i_segment + 1

        mesh.normals_split_custom_set(cn_list)

    """
    Basic paraxial lens analysis
    """
    #compute optical parameters
    y, u = 1, 0
    t_list = CT_list[1:]
    n_list = [1., ior1, ior2, ior3, 1.]
    n_elements = int(md)
    y1, u1 = paraxial.trace_lens(y,u, srad_list, t_list, n_list, n_elements)
    EFL = paraxial.calc_EFL(y, u1)
    self.flen_intern = EFL
            
    if display_edit:
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
    else:
        # do this explicitly in case edit was left toggled before   
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    
    return verts_outline
        

def add_rayfan(self, context):
    squarelens = self.squarelens
    ltype1 = self.ltype1
    ltype2 = self.ltype2
    ltype3 = self.ltype3
    ltype4 = self.ltype4
    RY1 = self.RY1
    RY2 = self.RY2
    RY3 = self.RY3
    RY4 = self.RY4
    RZ1 = self.RZ1
    RZ2 = self.RZ2
    RZ3 = self.RZ3
    RZ4 = self.RZ4
    surfrot1 = self.surfrot1
    surfrot2 = self.surfrot2
    surfrot3 = self.surfrot3
    surfrot4 = self.surfrot4
    N1 = self.num1
    N2 = self.num2
    lrad = self.lensradius
    flrad1 = self.flangerad1
    flrad2 = self.flangerad2
    flrad3 = self.flangerad3
    flrad4 = self.flangerad4
    CT1 = self.centerthickness1
    CT2 = self.centerthickness2
    CT3 = self.centerthickness3
    k1 = self.k1
    A1 = [a for a in self.A1]
    k2 = self.k2
    A2 = [a for a in self.A2]
    k3 = self.k3
    A3 = [a for a in self.A3]
    k4 = self.k4
    A4 = [a for a in self.A4]
    md = self.makedoublet
    dshape = self.dshape
    material_name1 = self.material_name1
    material_name2 = self.material_name2
    material_name3 = self.material_name3
    shade_smooth = self.shade_smooth
    smooth_type = self.smooth_type
    ior1 =  self.ior1
    ior2 =  self.ior2
    ior3 =  self.ior3
    display_edit = self.display_edit
    ghost_order = self.ghost_order

    ret_dummy = utils.check_surface(lrad, flrad1, RY1, k1, A1, RZ1, None, None, ltype1, squarelens)
    lrad1, hasfl1, flrad1, ltype1, surf_subtype_X1, surf_subtype_Y1, RY1, k1, A1, dSurfrot1 = ret_dummy
    ret_dummy = utils.check_surface(lrad, flrad2, RY2, k2, A2, RZ2, None, None, ltype2, squarelens)
    lrad2, hasfl2, flrad2, ltype2, surf_subtype_X2, surf_subtype_Y2, RY2, k2, A2, dSurfrot2 = ret_dummy
    ret_dummy = utils.check_surface(lrad, flrad3, RY3, k3, A3, RZ3, None, None, ltype3, squarelens)
    lrad3, hasfl3, flrad3, ltype3, surf_subtype_X3, surf_subtype_Y3, RY3, k3, A3, dSurfrot3 = ret_dummy
    ret_dummy = utils.check_surface(lrad, flrad4, RY4, k4, A4, RZ4, None, None, ltype4, squarelens)
    lrad4, hasfl4, flrad4, ltype4, surf_subtype_X4, surf_subtype_Y4, RY4, k4, A4, dSurfrot4 = ret_dummy
    
    surfrot1 = surfrot1 + dSurfrot1
    surfrot2 = surfrot2 + dSurfrot2
    surfrot3 = surfrot3 + dSurfrot3
    surfrot4 = surfrot4 + dSurfrot4


    surftype1 = surftype_Lens2Element(ltype1, surf_subtype_X1)
    surftype2 = surftype_Lens2Element(ltype2, surf_subtype_X2)
    surftype3 = surftype_Lens2Element(ltype3, surf_subtype_X3)
    surftype4 = surftype_Lens2Element(ltype4, surf_subtype_X4)

    # if not aspehric, unset k values
    # if not ltype1 == 'aspheric': k1 = 0
    # if not ltype2 == 'aspheric': k2 = 0
    # if not ltype3 == 'aspheric': k3 = 0
    # if not ltype4 == 'aspheric': k4 = 0
    ##################
    srad_list = [RY1, RY2, RY3, RY4]
    RZ_list = [RZ1, RZ2, RZ3, RZ4]
    lrad_list = [lrad1, lrad2, lrad3, lrad4]
    CT_list = [0, CT1, CT2, CT3, None] # Add None for handling in Element.add_surface()
    ltype_list = [ltype1, ltype2, ltype3, ltype4]
    lsubtype_list = [surf_subtype_X1, surf_subtype_X2, surf_subtype_X3, surf_subtype_X4]
    surftype_list = [surftype1, surftype2, surftype3, surftype4]
    surfrot_list = [surfrot1, surfrot2, surfrot3, surfrot4]
    k_list = [k1, k2, k3, k4]
    A_list = [A1, A2, A3, A4]
    n_list = [ior1, ior2, ior3]
    ##################
    n_surfaces = int(md) + 1
    y, u = 1, 0
    t_list = CT_list[1:-1] # "":-1" to exlcude sensor surface
    n_list2 = [1., ior1, ior2, ior3, 1.]
    n_elements = int(md)
    y1, u1 = paraxial.trace_lens(y,u, srad_list, t_list, n_list2, n_elements)
    BFL = paraxial.calc_BFL(y1, u1)
    ##################
   
    # create element and fill
    ele = Element()
    if squarelens: ele.outline_shape = 'square'
    for i in range(n_surfaces):
        #lrad1, hasfl1, flrad1, ssig1, surf_subtype1 = utils.check_surface(RY1, lrad, flrad1, k1, A1, ltype1)
        ele.add_surface(surftype_list[i], radius=srad_list[i], radius2=RZ_list[i], asph=[k_list[i]] + A_list[i],
                        rCA = lrad_list[i], CT=CT_list[i+1], surf_rotation=surfrot_list[i],
                        material = 'FIXVALUE_' + str(n_list2[i+1]))

   # create lens from element
    lens = Lenssystem()
    lens.elements = [[ele,0]]
    lens.num_optical_surfaces = n_surfaces
    if BFL >= 0:
        lens.detector['distance'] = BFL + self.zdet
    else:
        lens.detector['distance'] = 3.*sum(ele.data['CT']) + self.zdet
    if abs(lens.detector['distance']) > 1000*sum(ele.data['CT']): # protect from unrealistically large focal lengths, f/100 is already unusual
        lens.detector['distance'] = 1000*sum(ele.data['CT'])
    lens.detector['sizex'] = 10240 # default
    lens.detector['sizey'] = 10240 # default
    lens.detector['npixx'] = 10240 # default
    lens.detector['npixy'] = 10240 # default
    lens.detector['pixelpitch'] = 1 # default
    lens.build(0.586)
    # set up the rays
    initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist, self.fanangle1, self.fanangle2, self.fanangle3] 
    rays = rayfan.RayFan(self.fantype, initparams, store_history=True)
    # try to parse ghost_oder
    ghost_order = self.ghost_order
    nos = lens.num_optical_surfaces
    custom_surflist_valid = False
    if not ghost_order == '':
        try:
            surflist = [int(_) for _ in ghost_order.split(',')]
            if np.any(np.array(surflist) > nos) or np.any(np.array(surflist) < 1): 
                print("Warning: Custom ghost list invalid. Tracing default!")
            elif surflist[0] == surflist[1]:
                pass 
            else:
                s0 = min(surflist[:2])
                s1 = max(surflist[:2])
                seq1 = [_ for _ in range(1, s1 + 1)] 
                seq2 = [_ for _ in range(s0 + 1, s1)][::-1]
                seq3 = [_ for _ in range(s0, nos + 1)] 
                surflist = seq1 + seq2 + seq3
                custom_surflist_valid = True
        except:
            print("Warning: Custom ghost list invalid. Tracing default!")
    if not custom_surflist_valid:
        surflist = [i for i in range(1, nos + 1)] # standard surface list for non-ghost trace
    trace_detector = not self.tracetoscene
    rays = exec_trace(lens, rays, surflist, trace_detector=trace_detector)
    if trace_detector and self.autofocus:
        retval = rays.autofocus(EFL=self.flen_intern)
        if retval is not None:
            offset_min, P_new = retval
            rays.O_history[len(rays.O_history)-1] = P_new
        else:
            print('AUTOFOCUS FAILED')
        
    # if trace to scene, perform the last trace
    if self.tracetoscene:
        rays = trace_to_scene(context, rays)

    # create the lines
    verts = []
    edges = []
    faces = []
    ocr = False # self.onlycompleterays

    n_history = len(rays.O_history)
    n_rays = rays.O_history[0].shape[0]
    # create single numpy.array from history. Shape is (n_surf, n_rays, 3)
    O_hist = np.array([rays.O_history[i] for i in range(n_history)])
    idx_complete = ~np.isnan(O_hist[:,:,0].sum(axis=0))

    # reduce array down to fully valid rows
    if ocr:
        O_hist = O_hist[:, idx_complete, :]
        n_rays = O_hist.shape[1] # update this variable
    # change axis convention from raytrace to Blender orientation
    O_hist[:, :, [0, 1, 2]] = O_hist[:, :, [2, 0, 1]]
    O_hist[:, :, 0] *= -1
    # always add the verts for the first surface (i.e. ray fan input)
    o1 = np.copy(O_hist[0, :, :]) # initialize first surface
    verts = [list(o1[j, :]) for j in range(n_rays)]
    vidx1 = [j for j in range(n_rays)] # at each surface, the index in the verts list that each point has been assinged
    # loop over pairs of surfaces
    for i in range(1, n_history):
        offset1 = len(verts)
        # get new surface
        o2 = np.copy(O_hist[i, :, :])

        valid_both = ~np.isnan(o1[:, 0]) & ~np.isnan(o2[:, 0])
        n_valid = sum(valid_both)
        valid_any = np.any(valid_both)

        #add verts and lines
        if valid_any:
            verts_new = [list(o2[j, :]) for j, j_valid in enumerate(valid_both) if j_valid]
            verts = verts + verts_new
            vidx2 = np.full(n_rays, -1) # invalid points == -1
            vidx2[valid_both] = range(offset1, offset1 + n_valid)
            edges_new = [[vidx1[j], vidx2[j]] for j, j_valid in enumerate(valid_both) if j_valid]
            edges = edges + edges_new
        else:
            # can't continue tracing no valid ray
            break

        # switch for next loop iteration
        o2 = np.copy(o1)
        vidx1 = vidx2
    

    # calculate rms spot size
    if trace_detector:    
        # get points in detector plane
        P = np.array(rays.O_history[n_history - 1])
        # get mean point
        Pmean = np.nanmean(P, axis=0)
        # get differences
        Pdiff = P - Pmean
        r2 = np.einsum('ij,ij->i',Pdiff,Pdiff)
        # calculate rms
        rmsspotsize = np.sqrt(np.nanmean(r2))
    else:
        rmsspotsize = float('NaN')

    mesh = bpy.data.meshes.new(name = 'OC_Rayfan')
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)
    
    return rmsspotsize