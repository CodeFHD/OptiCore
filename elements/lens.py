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

class OBJECT_OT_add_lens(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "Lens"
    bl_options = {'REGISTER', 'UNDO'}
    
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
    ltype3 : EnumProperty(
           name="Surface 3 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 3",
           #options={'HIDDEN'},
           )
    ltype4 : EnumProperty(
           name="Surface 4 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 4",
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
           default = -24.,
           description="Radius of Curvature of Surface 2",
           unit = "LENGTH",
           )
    rad3 : FloatProperty(
           name="Surface 3 Radius",
           default = 24.,
           description="Radius of Curvature of Surface 3",
           unit = "LENGTH",
           )
    rad4 : FloatProperty(
           name="Surface 4 Radius",
           default = 24.,
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
           default = 3.,
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
           default = 1.,
           description="Center thickness of lens",
           unit = "LENGTH",
           )
    centerthickness2 : FloatProperty(
           name="Center Thickness 2",
           default = 1.,
           description="Center thickness of lens segment 2",
           unit = "LENGTH",
           )
    centerthickness3 : FloatProperty(
           name="Center Thickness 3",
           default = 1.,
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
           #default = (0.,0.,0.),
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
    material_name : StringProperty(
            name="Material",
            default="",
           )
    material_name2 : StringProperty(
            name="Material 2",
            default="",
           )
    material_name3 : StringProperty(
            name="Material 3",
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
            name="D-Shaped Lens",
            default=False,
           )
    optiverts : BoolProperty(
            name="Alternative radial vertex distribution",
            default=False,
           )
    debugmode : BoolProperty(
            name="Display Edit Mode",
            default=False,
           )
    ior1 : FloatProperty(
           name="IOR 1",
           default = 1.5,
           description="Index of Refraction of first lens section. Not transferred to material, only for paraxial infromations.",
           )
    ior2 : FloatProperty(
           name="IOR 2",
           default = 1.6,
           description="Index of Refraction of second lens section. Not transferred to material, only for paraxial infromations.",
           )
    ior3 : FloatProperty(
           name="IOR 3",
           default = 1.55,
           description="Index of Refraction of third lens section. Not transferred to material, only for paraxial infromations.",
           )
    def get_flen(self):
        return self.flen_intern
    flen_intern : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Paraxial EFL of the lens.",
           )
    flen : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Paraxial EFL of the lens.",
           get=get_flen,
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
    """
    fantype : EnumProperty(
           name="Ray Fan Type",
           items = {("f2d","2D",""),
                    ("f3d","3D",""),
                    ("f3dt","3D tris",""),
                    ("f3dr","3D random",""),
                    ("f2df","2D Finite",""),},
           default = "f2d",
           description="Ray Fan Type",
           #options={'HIDDEN'},
           )
    """
    fantype : EnumProperty(
           name="Ray Fan Type",
           items = {("2D","2D",""),
                    ("3D_square","3D square",""),
                    ("3D_tri","3D tris",""),
                    ("3D_random","3D random",""),
                    ("2D_finite","2D Finite",""),},
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
           min = 0.,
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

    def draw(self, context):
        md = self.makedoublet
        scene = context.scene
        #context.window_manager.invoke_props_dialog(self, width=500)
        row = self.layout
        #layout.scale_x = 2.0
        #row = layout.row(align=True)
        #row.scale_x = 2.0
        # Location
        col = row.column(align=True)
        col.label(text="Location")
        col.prop(self, 'location', text="")
        # Rotation
        col2 = col.column(align=True)
        col2.label(text="Rotation")
        col2.prop(self, 'rotation', text="")
        col3 = row.column(align=True)
        #col3.scale_x = 2.0
        col3.label(text="Lens Parameters")
        col3.prop(self, 'ltype1')
        col3.prop(self, 'rad1')
        if not self.ltype1=='aspheric': #TMP
            col3.prop(self, 'flangerad1')
        col3.prop(self, 'ltype2')
        col3.prop(self, 'rad2')
        if not self.ltype2=='aspheric': #TMP
            col3.prop(self, 'flangerad2')
        if md in ["2", "3"]:
            col3.prop(self, 'ltype3')
            col3.prop(self, 'rad3')
            if not self.ltype3=='aspheric': #TMP
                col3.prop(self, 'flangerad3')
        if md in ["3"]:
            col3.prop(self, 'ltype4')
            col3.prop(self, 'rad4')
            if not self.ltype4=='aspheric': #TMP
                col3.prop(self, 'flangerad4')
        col3.prop(self, 'lensradius')
        col3.prop(self, 'centerthickness1')
        if md in ["2", "3"]:
            col3.prop(self, 'centerthickness2')
        if md in ["3"]:
            col3.prop(self, 'centerthickness3')
        col3.prop(self, 'num1')
        col3.prop(self, 'num2')
        if self.ltype1=='aspheric':
            col3.prop(self, 'k1')
            col3.prop(self, 'A1')
        if self.ltype2=='aspheric':
            col3.prop(self, 'k2')
            col3.prop(self, 'A2')
        if md in ["2", "3"]:
            if self.ltype3=='aspheric':
                col3.prop(self, 'k3')
                col3.prop(self, 'A3')
        if md in ["3"]:
            if self.ltype4=='aspheric':
                col3.prop(self, 'k4')
                col3.prop(self, 'A4')
        # col3.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        # if md in ["2", "3"]:
        #     col3.prop_search(self, "material_name2", bpy.data, "materials", icon="NONE")
        # if md in ["3"]:
        #     col3.prop_search(self, "material_name3", bpy.data, "materials", icon="NONE")
        col3.prop(self, 'makedoublet')
        col3.prop(self, 'shade_smooth')
        # if self.shade_smooth:
        #     col3.prop(self, 'smooth_type')
        col3.prop(self, 'dshape')
        #col3.prop(self, 'optiverts')
        col3.prop(self, 'debugmode')
        col4 = row.column(align=True)
        #col4.scale_x = 2.0
        col4.label(text="Optical Parameters")
        col4.prop(self, 'ior1')
        if md in ["2", "3"]:
            col4.prop(self, 'ior2')
        if md in ["3"]:
            col4.prop(self, 'ior3')
        col4.prop(self, 'flen')
        col4.prop(self, 'addrayfan')
        if self.addrayfan:
            col4.prop(self, 'zdet')
            col4.prop(self, 'nrays')
            col4.prop(self, 'fantype')
            col4.prop(self, 'fandist')
            col4.prop(self, 'fanangle')
            col4.prop(self, 'fandiam')
            col4.prop(self, 'tracetoscene')

    def execute(self, context):
        add_lens(self, context)
        if self.addrayfan:
        #     utils.trace_rays(self, context)
            add_rayfan(self, context)
        return {'FINISHED'}
    
def get_default_paramdict():
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
    paramdict['material_name'] = ""
    paramdict['material_name2'] = ""
    paramdict['material_name3'] = ""
    paramdict['shade_smooth'] = True
    paramdict['smooth_type'] = True
    paramdict['ior1'] = 1.5
    paramdict['ior2'] = 1.6
    paramdict['ior3'] = 1.55
    paramdict['debugmode'] = False
    return paramdict

def add_lens(self, context, paramdict=None):
    edges = []
    
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
        material_name = self.material_name
        material_name2 = self.material_name2
        material_name3 = self.material_name3
        shade_smooth = self.shade_smooth
        smooth_type = self.smooth_type
        ior1 =  self.ior1
        ior2 =  self.ior2
        ior3 =  self.ior3
        debugmode = self.debugmode
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
        material_name = paramdict['material_name']
        material_name2 = paramdict['material_name2']
        material_name3 = paramdict['material_name3']
        shade_smooth = paramdict['shade_smooth']
        smooth_type = paramdict['smooth_type']
        ior1 =  paramdict['ior1']
        ior2 =  paramdict['ior2']
        ior3 =  paramdict['ior3']
        debugmode = paramdict['debugmode']
    
    num_sides = int(md)       
    num_surfaces = int(md) + 1

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
    
    #check surface radii for consistency
    ##check radius overflow
    if not utils.check_surface(np.abs(srad1), lrad1):
        srad1 = 0
        lrad1 = lrad
        hasfl1 = 0
    if not utils.check_surface(np.abs(srad2), lrad2):
        srad2 = 0
        lrad2 = lrad
        hasfl2 = 0
    if md in ["2", "3"]:
        if not utils.check_surface(np.abs(srad3), lrad3):
            srad3 = 0
            lrad3 = lrad
            hasfl3 = 0
    if md in ["3"]:
        if not utils.check_surface(np.abs(srad4), lrad4):
            srad4 = 0
            lrad4 = lrad
            hasfl4 = 0
    hasfl = [hasfl1, hasfl2, hasfl3, hasfl4]
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
        
    # collect into lists for looping
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

    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces, normals = sfc.add_flat_surface(lrad1,N1,N2,1, dshape=dshape)
    elif ltype1 == 'spherical':
        verts, faces, normals, N1, N2 = sfc.add_spherical_surface(srad1, lrad1, N1, N2,1, dshape=dshape, optiverts=optiverts, lrad_ext=lrad)
    elif ltype1 == 'aspheric':
        verts, faces, normals = sfc.add_aspheric_surface(srad1, k1, A1, lrad1, N1, N2,1, dshape=dshape)
    
    nVerts1 = len(verts)
    nFacs1 = len(faces)
    nVerts_at_segment.append(len(verts))
    
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
        # they are named with index 0 to avoid confusion with globals without number where exits, e.g. lrad
        srad0 = srad_list[i]
        lrad0 = lrad_list[i]
        CT_sum = CT_sum + CT_list[i]
        ltype0 = ltype_list[i]
        k0 = k_list[i]
        A0 = A_list[i]
        
        # get the new surface and append
        if srad0 == 0:#flat surface case
            dvert, dfac, dnormals = sfc.add_flat_surface(lrad0,N1,N2,+1,CT_sum,nVerts=nVerts_tot, dshape=dshape)
        elif ltype0 == 'spherical':
            dvert, dfac, dnormals, N12, N22 = sfc.add_spherical_surface(srad0, lrad0, N1, N2, +1, CT_sum, nVerts=nVerts_tot, dshape=dshape, optiverts=optiverts, lrad_ext=lrad)
        elif ltype0 == 'aspheric':
            dvert, dfac, dnormals = sfc.add_aspheric_surface(srad0, k0, A0, lrad0, N1, N2, +1, CT_sum, nVerts=nVerts_tot, dshape=dshape)
        # flip normals for last surface
        if i==int(md):
            dfac = [df[::-1] for df in dfac]
            dnormals = [(-n[0], -n[1], -n[2]) for n in dnormals]

        verts = verts + dvert
        faces = faces + dfac
        normals = normals + dnormals
        
        nVerts_this = len(dvert)
        nVerts_at_segment.append(len(verts))
        
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
                #ptss1 = [nVerts_tot - (N2-1), nVerts_tot]
            else:
                p0 = nVerts_tot - nVerts_prev 
                ptss1 = [p0 + (x+1)*N2 for x in range(N1)][::-1] + [p0] + [p0 + 1 + x*N2 for x in range(N1)]
            if srad0==0:
                ptss2 = [len(verts)-N2, len(verts)-1]
            else:
                p0 = len(verts) - nVerts_this
                ptss2 = [p0 + 1 + x*N2 for x in range(N1)][::-1] + [p0] + [p0 + (x+1)*N2 for x in range(N1)]
                #ptss2 = [nVerts_tot + x*N2 for x in range(N1)] + [len(verts)-1] + [nVerts_tot + (x+1)*N2-1 for x in range(N1)[::-1]]
            faces.append(ptss1 + ptss2)
        
        # update needed carry-over-parameters after this surface
        nVerts_prev = len(dvert)
        #if i == 1: nVerts_prev = N2
        nVerts_tot = len(verts) # the (total) length of verts so far
        nVerts_at_segment.append(len(verts))
        srad_prev = 1.*srad0
    
    #create mesh from verts and faces
    del dvert
    del dfac
    mesh = bpy.data.meshes.new(name="Lens")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    """
    #assign material(s)
    if material_name in bpy.data.materials:
        mat = bpy.data.materials[material_name]
        obj.data.materials.append(mat)
    if md:
        cond2 = material_name2 in bpy.data.materials
        cond3 = material_name3 in bpy.data.materials
    if md and cond2 and cond3:
        if material_name2 in bpy.data.materials:
            mat2 = bpy.data.materials[material_name2]
            obj.data.materials.append(mat2)
        if material_name3 in bpy.data.materials:
            mat3 = bpy.data.materials[material_name3]
            obj.data.materials.append(mat3)
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [False, False, True]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        #rear surface (mat3)
        for i in range(nFacSide):
            mesh.polygons[i + nFacs3].select=True
        for i in range(ndFac3):
            mesh.polygons[i + nFacs2+ nFacSide].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        obj.active_material_index = 2
        bpy.ops.object.material_slot_assign()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        #middle surface (mat2)
        for i in range(ndFac2):
            mesh.polygons[i + nFacs].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        obj.active_material_index = 1
        bpy.ops.object.material_slot_assign()
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    """
    
    # apply smooth shading
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

    #compute optical parameters
    y, u = 1, 0
    t_list = CT_list[1:]
    n_list = [1., ior1, ior2, ior3, 1.]
    n_elements = int(md)
    y1, u1 = paraxial.trace_lens(y,u, srad_list, t_list, n_list, n_elements)
    EFL = paraxial.calc_EFL(y, u1)
    self.flen_intern = EFL
    # if not md:
    #     self.flen_intern = utils.lens_math.f_lensmaker(srad1, -srad2, ior1, CT1)
            
    if debugmode:
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        

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
    material_name = self.material_name
    material_name2 = self.material_name2
    material_name3 = self.material_name3
    shade_smooth = self.shade_smooth
    smooth_type = self.smooth_type
    ior1 =  self.ior1
    ior2 =  self.ior2
    ior3 =  self.ior3
    debugmode = self.debugmode
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
    ##################
    srad_list = [srad1, srad2, srad3, srad4]
    lrad_list = [lrad1, lrad2, lrad3, lrad4]
    CT_list = [0, CT1, CT2, CT3]
    ltype_list = [ltype1, ltype2, ltype3, ltype4]
    k_list = [k1, k2, k3, k4]
    #A_list = [A1, A2, A3, A4]
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
    ele.data['type'] = ['STANDARD' for i in range(n_surfaces)]
    ele.data['radius'] = srad_list[:n_surfaces]
    ele.data['asph'] = [None for i in range(n_surfaces)]
    ele.data['CT'] = CT_list[1:n_surfaces]
    # ele.data['lrad'].append(surf_infos[idx]['lrad'])
    ele.data['rCA'] = lrad_list[:n_surfaces]
    ele.data['decenter'] = [None for i in range(n_surfaces)]
    ele.data['tilt'] = [None for i in range(n_surfaces)]
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
        initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist, -1.*rayfany]
    else:
        initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist]    
    rays = rayfan.RayFan(self.fantype, initparams, store_history=True)
    surflist = [i for i in range(1, lens.num_optical_surfaces + 1)] # standard surface list for non-ghost trace
    trace_detecor = not self.tracetoscene
    rays = exec_trace(lens, rays, surflist, trace_detector=trace_detecor)
        
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
            
    mesh = bpy.data.meshes.new(name="Rayfan")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)