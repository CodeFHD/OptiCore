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

import bpy
import numpy as np

from mathutils import Vector

from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils
from ..utils.paraxial import paraxial
from ..utils.raytrace.element import Element
from ..utils.raytrace.lenssystem import Lenssystem
from ..utils.raytrace.trace_sequential import exec_trace, trace_to_scene
from ..utils.raytrace import rayfan

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
        items = {("default","Most Important",""),
                 ("optical","Optical Design",""),
                 ("geometry","Model Geometry",""),
                 ("rayfan","Ray Tracing",""),
                 ("paraxial","Paraxial Optics",""),},
        default = "default",
        description="Select which settings are displayed in this window.",
           )
    makedoublet : EnumProperty(
        name="Lens configuration",
        items = {("1", "Singlet", ""),
                 ("2", "Doublet", ""),
                 ("3", "Triplet", "")},
        default="1",
        )
    ltype1 : EnumProperty(
           name="Surface 1 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "aspheric",
           description="Shape of Surface 1",
           #options={'HIDDEN'},
           )
    ltype2 : EnumProperty(
           name="Surface 2 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "aspheric",
           description="Shape of Surface 2",
           #options={'HIDDEN'},
           )
    ltype3 : EnumProperty(
           name="Surface 3 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "aspheric",
           description="Shape of Surface 3",
           #options={'HIDDEN'},
           )
    ltype4 : EnumProperty(
           name="Surface 4 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "aspheric",
           description="Shape of Surface 4",
           #options={'HIDDEN'},
           )
    rad1 : FloatProperty(
           name="Surface 1 Radius",
           default = 61.47,
           description="Radius of Curvature of Surface 1",
           unit = "LENGTH",
           )
    rad2 : FloatProperty(
           name="Surface 2 Radius",
           default = -44.64,
           description="Radius of Curvature of Surface 2",
           unit = "LENGTH",
           )
    rad3 : FloatProperty(
           name="Surface 3 Radius",
           default = -129.94,
           description="Radius of Curvature of Surface 3",
           unit = "LENGTH",
           )
    rad4 : FloatProperty(
           name="Surface 4 Radius",
           default = -126.,
           description="Radius of Curvature of Surface 4",
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
           default = 12.5,
           description="Lens outer radius",
           unit = "LENGTH",
           )
    flangerad1 : FloatProperty(
           name="Flange Width Surface 1",
           default = 0.,
           description="Flange Width Surface 1",
           unit = "LENGTH",
           min = 0,
           )
    flangerad2 : FloatProperty(
           name="Flange Width Surface 2",
           default = 0.,
           description="Flange Width Surface 2",
           unit = "LENGTH",
           min = 0,
           )
    flangerad3 : FloatProperty(
           name="Flange Width Surface 3",
           default = 0.,
           description="Flange Width Surface 3",
           unit = "LENGTH",
           min = 0,
           )
    flangerad4 : FloatProperty(
           name="Flange Width Surface 4",
           default = 0.,
           description="Flange Width Surface 4",
           unit = "LENGTH",
           min = 0,
           )
    centerthickness1 : FloatProperty(
           name="Center Thickness",
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
           name="k",
           default = 0.,
           description="Aspheric conical constant",
           )
    A1 : FloatVectorProperty(
           name="A",
           # default = (1e-5,0.,0.), # for testing
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    k2 : FloatProperty(
           name="k2",
           default = 0.,
           description="Aspheric conical constant",
           )
    A2 : FloatVectorProperty(
           name="A2",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    k3 : FloatProperty(
           name="k3",
           default = 0.,
           description="Aspheric conical constant",
           )
    A3 : FloatVectorProperty(
           name="A3",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    k4 : FloatProperty(
           name="k4",
           default = 0.,
           description="Aspheric conical constant",
           )
    A4 : FloatVectorProperty(
           name="A4",
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
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
    optiverts : BoolProperty(
            name="Alternative radial vertex distribution",
            default=False,
           )
    display_edit : BoolProperty(
            name="Display Edit Mode",
            default=False,
           )
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
           items = {("2D","2D",""),
                    ("3D_square","3D square",""),
                    ("3D_tri","3D tris",""),
                    ("3D_random","3D random",""),
                    ("2D_finite","2D Finite",""),
                    ("3D_rings", "3D rings", ""),},
           default = "2D",
           description="Ray Fan Type",
           #options={'HIDDEN'},
           )
    fandist : FloatProperty(
           name="Ray Fan Origin",
           default = 20.,
           description="Distance where ray fan originates.",
           )
    fanangle : FloatProperty(
           name="Ray Fan Angle",
           default = 0.,
           description="Angle of Ray Fan.",
           min = -90.,
           max = 90.,
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
        items = {("r1","Surface 1 Radius",""),
                 ("r2","Surface 2 Radius",""),
                 ("r3","Surface 3 Radius",""),
                 ("r4","Surface 4 Radius",""),
                 ("n1","Glass 1 index",""),
                 ("n2","Glass 2 index",""),
                 ("n3","Glass 3 index",""),},
        default = "r1",
        description="Select which parameter to automatically optimize.",
           )
    flen_goal : FloatProperty(
           name="Focal Length target",
           default = 100.,
           description="Paraxial EFL of the lens used in autosolver.",
           min=1.,
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
            col.prop(self, 'rad1')
            col.prop(self, 'rad2')
            col.prop(self, 'rad3')
            col.prop(self, 'rad4')
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
            col.prop(self, 'centerthickness1')
            col.prop(self, 'centerthickness2')
            col.prop(self, 'centerthickness3')
            # col.prop(self, 'ltype1')
            col.prop(self, 'rad1')
            col.prop(self, 'k1')
            col.prop(self, 'A1')
            col.prop(self, 'flangerad1')
            # col.prop(self, 'ltype2')
            col.prop(self, 'rad2')
            col.prop(self, 'k2')
            col.prop(self, 'A2')
            col.prop(self, 'flangerad2')
            # col.prop(self, 'ltype3')
            col.prop(self, 'rad3')
            col.prop(self, 'k3')
            col.prop(self, 'A3')
            col.prop(self, 'flangerad3')
            # col.prop(self, 'ltype4')
            col.prop(self, 'rad4')
            col.prop(self, 'k4')
            col.prop(self, 'A4')
            col.prop(self, 'flangerad4')
        if disp == 'geometry':
            # Location
            col.label(text="Location")
            col.prop(self, 'location', text="")
            # Rotation
            col.label(text="Rotation")
            col.prop(self, 'rotation', text="")
            col.label(text="Modelling Parameters")
            col.prop(self, 'num1')
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
            col.prop(self, 'fanangle')
            col.prop(self, 'fandiam')
            col.prop(self, 'tracetoscene')
            col.prop(self, 'autofocus')
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
    paramdict['ltype1'] = "spherical"
    paramdict['ltype2'] = "spherical"
    paramdict['ltype3'] = "spherical"
    paramdict['ltype4'] = "spherical"
    paramdict['rad1'] = 12.
    paramdict['rad2'] = 24.
    paramdict['rad3'] = 24.
    paramdict['rad4'] = 24.
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
    paramdict['optiverts'] = False
    paramdict['material_name1'] = ""
    paramdict['material_name2'] = ""
    paramdict['material_name3'] = ""
    paramdict['material_name4'] = ""
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
        ltype1 = self.ltype1
        ltype2 = self.ltype2
        ltype3 = self.ltype3
        ltype4 = self.ltype4
        srad1 = self.rad1
        srad2 = self.rad2
        srad3 = self.rad3
        srad4 = self.rad4
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
        optiverts = self.optiverts
        material_name1 = self.material_name1
        material_name2 = self.material_name2
        material_name3 = self.material_name3
        material_name4 = self.material_name4
        shade_smooth = self.shade_smooth
        smooth_type = self.smooth_type
        ior1 =  self.ior1
        ior2 =  self.ior2
        ior3 =  self.ior3
        display_edit = self.display_edit
    else:
        ltype1 = paramdict['ltype1']
        ltype2 = paramdict['ltype2']
        ltype3 = paramdict['ltype3']
        ltype4 = paramdict['ltype4']
        srad1 = paramdict['rad1']
        srad2 = paramdict['rad2']
        srad3 = paramdict['rad3']
        srad4 = paramdict['rad4']
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
        optiverts = paramdict['optiverts']
        material_name1 = paramdict['material_name1']
        material_name2 = paramdict['material_name2']
        material_name3 = paramdict['material_name3']
        material_name4 = paramdict['material_name4']
        shade_smooth = paramdict['shade_smooth']
        smooth_type = paramdict['smooth_type']
        ior1 =  paramdict['ior1']
        ior2 =  paramdict['ior2']
        ior3 =  paramdict['ior3']
        display_edit = paramdict['display_edit']

    """
    Verification of the geometry
    """

    lrad1, hasfl1, flrad1, ssig1 = utils.check_surface(srad1, lrad, flrad1, k1, A1)
    lrad2, hasfl2, flrad2, ssig2 = utils.check_surface(srad2, lrad, flrad2, k2, A2)
    lrad3, hasfl3, flrad3, ssig3 = utils.check_surface(srad3, lrad, flrad3, k3, A3)
    lrad4, hasfl4, flrad4, ssig4 = utils.check_surface(srad4, lrad, flrad4, k4, A4)
   
    ##check center thickness
    """        
    lsurf1, lsurf2 = 0, 0
    if not srad1 == 0:
        lsurf1 = srad1-ssig1*np.sqrt(srad1**2-lrad1**2)
    if not srad2 == 0:
        lsurf2 = srad2-ssig2*np.sqrt(srad2**2-lrad2**2)
    if (lsurf1 + lsurf2) > CT1:
        CT1 = lsurf1 + lsurf2
    if md in ["2", "3"]:
        lsurf3 = 0
        if not srad3 == 0:
            lsurf3 = srad3-ssig3*np.sqrt(srad3**2-lrad3**2)
        if (-lsurf2 + lsurf3) > CT2: CT2 = -lsurf2 + lsurf3
    if md in ["3"]:
        lsurf4 = 0
        if not srad4 == 0:
            lsurf4 = srad4-ssig4*np.sqrt(srad4**2-lrad4**2)
        if (-lsurf3 + lsurf4) > CT3: CT3 = -lsurf3 + lsurf4
    """
        
    """
    prepare some variables
    """
    num_sides = int(md)       
    num_surfaces = int(md) + 1
    
    # collect into lists for looping
    hasfl = [hasfl1, hasfl2, hasfl3, hasfl4]
    srad_list = [srad1, srad2, srad3, srad4]
    lrad_list = [lrad1, lrad2, lrad3, lrad4]
    CT_list = [0, CT1, CT2, CT3]
    ltype_list = [ltype1, ltype2, ltype3, ltype4]
    k_list = [k1, k2, k3, k4]
    A_list = [A1, A2, A3, A4]
    
    # for keeping up with counts
    # named segment because sufacesd and sides alternate S1, S2, E1, S3, E2, ...
    # initialize with zero in order to facilitate taking differences starting from surface 0
    nVerts_at_segment = [0]
    nLoops_at_segment = [0]
    nFaces_at_segment = [0]
    
    """
    Generating the vertices and faces
    """
    
    # keeping the separate generator algorithms for consistency testing and efficiency
    if k1 == 0 and np.all(np.array(A1) == 0):
        ltype1 = 'spherical'

    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces, normals = sfc.add_flat_surface(lrad1,N1,N2,1, dshape=dshape)
    elif ltype1 == 'spherical':
        verts, faces, normals, N1, N2 = sfc.add_spherical_surface(srad1, lrad1, N1, N2,1, dshape=dshape, optiverts=optiverts, lrad_ext=lrad)
    elif ltype1 == 'aspheric':
        verts, faces, normals = sfc.add_aspheric_surface(srad1, k1, A1, lrad1, N1, N2,1, dshape=dshape, lrad_ext=lrad)
    
    nVerts1 = len(verts)
    nFacs1 = len(faces)
    nVerts_at_segment.append(len(verts))
    nFaces_at_segment.append(len(faces))
    
    # all sides are equal so this is generic and gets added below
    normalsside = sfc.get_ringnormals(N2, dshape=dshape)
    nFacSide = N2
    
    # add further surfaces
    CT_sum = 0
    srad_prev = srad1
    nVerts_tot = nVerts1
    nVerts_prev = nVerts1 # N2
    for i in range(1, int(md)+1):
        # get parameters for this surface
        # they are named with index 0 to avoid confusion with globals without number where that exists, e.g. lrad
        srad0 = srad_list[i]
        lrad0 = lrad_list[i]
        CT_sum = CT_sum + CT_list[i]
        ltype0 = ltype_list[i]
        k0 = k_list[i]
        A0 = A_list[i]
        if k0 == 0 and np.all(np.array(A0) == 0):
            ltype0 = 'spherical'
        
        # get the new surface and append
        if srad0 == 0:#flat surface case
            dvert, dfac, dnormals = sfc.add_flat_surface(lrad0,N1,N2,+1,CT_sum,nVerts=nVerts_tot, dshape=dshape)
        elif ltype0 == 'spherical':
            dvert, dfac, dnormals, N12, N22 = sfc.add_spherical_surface(srad0, lrad0, N1, N2, +1, CT_sum, nVerts=nVerts_tot, dshape=dshape, optiverts=optiverts, lrad_ext=lrad)
        elif ltype0 == 'aspheric':
            dvert, dfac, dnormals = sfc.add_aspheric_surface(srad0, k0, A0, lrad0, N1, N2, +1, CT_sum, nVerts=nVerts_tot, dshape=dshape, lrad_ext=lrad)
        # flip normals for last surface
        if i==int(md):
            dfac = [df[::-1] for df in dfac]
            dnormals = [(-n[0], -n[1], -n[2]) for n in dnormals]

        verts = verts + dvert
        faces = faces + dfac
        normals = normals + dnormals
        
        nVerts_this = len(dvert)
        nVerts_at_segment.append(len(verts))
        nFaces_at_segment.append(len(faces))
        
        # add the side
        for j in range(N2-dshape):
            fi1 = len(verts) - N2 + (j+1)%(N2)
            fi2 = len(verts) - N2 + j
            fi3 = fi2 - nVerts_this
            fi4 = fi1 - nVerts_this
            faces.append([fi1,fi2,fi3,fi4])
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
            faces.append(ptss1 + ptss2)
        
        # update needed carry-over-parameters after this surface
        nVerts_prev = len(dvert)
        nVerts_tot = len(verts) # the (total) length of verts so far
        nVerts_at_segment.append(len(verts))
        nFaces_at_segment.append(len(faces))
        srad_prev = 1.*srad0
    
    #create mesh from verts and faces
    del dvert
    del dfac
    edges = [] # edges are not explicitly created here so we pass an empty list
    mesh = bpy.data.meshes.new(name="Lens")
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
    hasmat_any = hasmat1 or hasmat2 or hasmat3 or hasmat4
    
    # initially change to edit mode
    if hasmat_any:
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [False, False, True] # face select mode
    
    # apend the materials to the objects material list
    dummy = 0 # keep track of the material index in the list if there is a gap
    if hasmat1:
        mat1 = bpy.data.materials[material_name1]
        obj.data.materials.append(mat1)
        mat_idx_1 = dummy
        dummy = dummy + 1
    if hasmat2:
        mat2 = bpy.data.materials[material_name2]
        obj.data.materials.append(mat2)
        mat_idx_2 = dummy
        dummy = dummy + 1
    if hasmat3:
        mat3 = bpy.data.materials[material_name3]
        obj.data.materials.append(mat3)
        mat_idx_3 = dummy
        dummy = dummy + 1
    if hasmat4:
        mat4 = bpy.data.materials[material_name4]
        obj.data.materials.append(mat4)
        mat_idx_4 = dummy
        
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
        obj.active_material_index = mat_idx_1
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
        obj.active_material_index = mat_idx_2
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
        obj.active_material_index = mat_idx_3
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
        obj.active_material_index = mat_idx_4
        bpy.ops.object.material_slot_assign()
    
    # side
    # TODO implement if requested
    
    # d-shaoe
    # TODO implement if requested
    
    """
    Creating split edge normals for smooth shading
    """
    
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    if shade_smooth:
        bpy.ops.object.shade_smooth()
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
        bpy.ops.mesh.customdata_custom_splitnormals_add()
        cn_list = []
        i_segment = 1 # start at 1 because of 0 offset
        n_loops = 0
        n_loops_prev = 0
        for i_surf in range(num_surfaces):
            # surfaces
            srad0 = srad_list[i_surf]
            if srad0 == 0:
                nloops1 = N2 # when radius == 0, the face is flat circle, only one loop line going around
                nloops1_fl = N2 # flange has no impact for a flat face
            else:
                nloops1 = (N2-dshape)*(3 + 4*(N1-1)) # center of lens has triangular faces, hence the 3. Rest quads, hence the 4. N-1 because of the innter triangle section.
                nloops1_fl = nloops1 - hasfl[i_surf]*4*(N2-dshape) # 2*hasfl[i_surf]*(3 + 4*(N1-1)+1)
            for i in range(nloops1_fl):
                vi = mesh.loops[i + n_loops].vertex_index
                cn_list.append(normals[vi])
            if hasfl[i_surf]:
                for i in range(4*(N2-dshape)): # (2*(3 + 4*(N1-1))+2):
                    vi = mesh.loops[nloops1_fl + n_loops + i].vertex_index
                    cn_list.append(normals[nVerts_at_segment[i_segment]-1]) # all annulus vertices have the same pointing
            i_segment = i_segment + 1
            n_loops_prev = n_loops
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
            n_loops_prev = n_loops
            n_loops = n_loops + n_loop_side
        
            # D-shape
            if dshape:
                if srad_list[i_surf-1] == 0:
                    n_side1 = 1
                else:
                    n_side1 = 2*N1
                if srad_list[i_surf]== 0:
                    n_side2 = 1
                else:
                    n_side2 = 2*N1
                for i in range(n_side1 + n_side2 + 2):
                    cn_list.append((0,-1,0))
                n_loops_prev = n_loops
                n_loops = n_loops + n_side1 + n_side2 + 2
            
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
        

def add_rayfan(self, context):
    ##################
    ltype1 = self.ltype1
    ltype2 = self.ltype2
    ltype3 = self.ltype3
    ltype4 = self.ltype4
    srad1 = self.rad1
    srad2 = self.rad2
    srad3 = self.rad3
    srad4 = self.rad4
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
    optiverts = self.optiverts
    material_name1 = self.material_name1
    material_name2 = self.material_name2
    material_name3 = self.material_name3
    shade_smooth = self.shade_smooth
    smooth_type = self.smooth_type
    ior1 =  self.ior1
    ior2 =  self.ior2
    ior3 =  self.ior3
    display_edit = self.display_edit
    ##################
    hasfl1 = flrad1 > 0.01*lrad
    hasfl2 = flrad2 > 0.01*lrad
    if ltype1 == 'aspheric':
        hasfl1 = 0 #TMP
    if ltype2 == 'aspheric':
        hasfl2 = 0 #TMP
    ssig1 = 1
    if srad1 < 0:
        ssig1 = -1
    ssig2 = 1
    if srad2 < 0:
        ssig2 = -1

    if md in ["2", "3"]:
        hasfl3 = flrad3 > 0.01*lrad
        if ltype3 == 'aspheric':
            hasfl3 = 0 #TMP
        ssig3 = 1
        if srad3 < 0:
            ssig3 = -1
    else:
        hasfl3 = 0
    if md in ["3"]:
        hasfl4 = flrad4 > 0.01*lrad
        if ltype4 == 'aspheric':
            hasfl4 = 0 #TMP
        ssig4 = 1
        if srad4 < 0:
            ssig4 = -1
    else:
        hasfl4 = 0
    #TODO: This may need to move down
    if flrad1 > 0.99*lrad: flrad1 = 0.99*lrad
    if flrad2 > 0.99*lrad: flrad2 = 0.99*lrad
    lrad1 = lrad - flrad1 if hasfl1 else lrad
    lrad2 = lrad - flrad2 if hasfl2 else lrad
    if md in ["2", "3"]:
        if flrad3 > 0.99*lrad: flrad3 = 0.99*lrad
        lrad3 = lrad - flrad3 if hasfl3 else lrad
    else:
        lrad3 = 0
    if md in ["3"]:
        if flrad4 > 0.99*lrad: flrad4 = 0.99*lrad
        lrad4 = lrad - flrad4 if hasfl4 else lrad
    else:
        lrad4 = 0
    # if not aspehric, unset k values
    if not ltype1 == 'aspheric': k1 = 0
    if not ltype2 == 'aspheric': k2 = 0
    if not ltype3 == 'aspheric': k3 = 0
    if not ltype4 == 'aspheric': k4 = 0
    ##################
    srad_list = [srad1, srad2, srad3, srad4]
    lrad_list = [lrad1, lrad2, lrad3, lrad4]
    CT_list = [0, CT1, CT2, CT3]
    ltype_list = [ltype1, ltype2, ltype3, ltype4]
    k_list = [k1, k2, k3, k4]
    A_list = [A1, A2, A3, A4]
    n_list = [ior1, ior2, ior3]
    ##################
    n_surfaces = int(md) + 1
    y, u = 1, 0
    t_list = CT_list[1:]
    n_list2 = [1., ior1, ior2, ior3, 1.]
    n_elements = int(md)
    y1, u1 = paraxial.trace_lens(y,u, srad_list, t_list, n_list2, n_elements)
    BFL = paraxial.calc_BFL(y1, u1)
    ##################
   
    # create element and fill
    ele = Element()
    # lists with length n_surfaces
    ele.data['type'] = ['STANDARD' for i in range(n_surfaces)]
    ele.data['radius'] = srad_list[:n_surfaces]
    ele.data['asph'] = [[k_list[i]] + A_list[i] for i in range(n_surfaces)]
    ele.data['coating'] = [None for i in range(n_surfaces)]
    # ele.data['lrad'].append(surf_infos[idx]['lrad'])
    ele.data['rCA'] = lrad_list[:n_surfaces]
    ele.data['surf_decenter'] = [None for i in range(n_surfaces)]
    ele.data['surf_tilt'] = [None for i in range(n_surfaces)]
    # lists with length n_surfaces - 1
    ele.data['CT'] = CT_list[1:n_surfaces]
    ele.data['material'] = ['FIXVALUE_' + str(n_list[i]) for i in range(n_surfaces-1)]
   
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
    lens.detector['pixelpitch'] = 0.01 # default
    lens.build(0.586)
    # set up the rays
    if self.fantype in ["2D_finite"]:
        rayfany = self.fandist*np.tan(self.fanangle*np.pi/180)
        initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist, 1.*rayfany]
    else:# self.fantype in ["2D"]:
        initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist, self.fanangle*np.pi/180] 
    #else:
    #    initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist]    
    rays = rayfan.RayFan(self.fantype, initparams, store_history=True)
    surflist = [i for i in range(1, lens.num_optical_surfaces + 1)] # standard surface list for non-ghost trace
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
    nVerts = 0
    for j in range(len(rays.O_history)-1):
        for i in range(rays.O_history[0].shape[0]):
            # change axis convention from raytrace to Blender orientation
            o = np.array(rays.O_history[j][i])
            o[[0,1,2]] = o[[2,0,1]]
            o[0] = -1*o[0]
            verts.append(Vector(o))
            o = np.array(rays.O_history[j+1][i])
            o[[0,1,2]] = o[[2,0,1]]
            o[0] = -1*o[0]
            verts.append(Vector(o))
            edges.append([nVerts + 2*i, nVerts + 2*i + 1])
        nVerts = len(verts)
        
    # calculate rms spot size
    if trace_detector:    
        # get points in detector plane
        P = np.array(rays.O_history[len(rays.O_history)-1])
        # get mean point
        Pmean = np.nanmean(P, axis=0)
        # get differences
        Pdiff = P - Pmean
        r2 = np.einsum('ij,ij->i',Pdiff,Pdiff)
        # calculate rms
        rmsspotsize = np.sqrt(np.nanmean(r2))
    else:
        rmsspotsize = float('NaN')
            
    mesh = bpy.data.meshes.new(name="Rayfan")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)
    
    return rmsspotsize