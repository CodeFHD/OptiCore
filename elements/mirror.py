import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils

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
        #scene = context.scene
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

        add_mirror(self, context)

        return {'FINISHED'}

def add_mirror(self, context):
    edges = []
    
    srad = self.rad
    N1 = self.num1
    N2 = self.num2
    mrad = self.mirrorradius
    CT = self.centerthickness
    surftype = self.mtype
    theta = self.theta*180/np.pi
    opos = self.opos
    hrad=0
    if self.cent_hole:
        hrad = self.hole_rad
        if hrad > 0.99*mrad:
            hrad = 0.99*mrad
    
    #check surface radius for consistency
    if surftype == 'spherical':
        if not utils.check_surface(np.abs(srad), mrad): srad=0

    #compute mirror surface
    if srad == 0: #flat surface case
        verts, faces, splitverts = sfc.add_flat_surface(mrad,N1,N2,hole=self.cent_hole,hrad=hrad)
        yadd = 0
        xOA = 0
    else:
        if surftype == 'spherical':
            verts, faces, splitverts = sfc.add_spherical_surface(-srad, mrad, N1, N2,hole=self.cent_hole,hrad=hrad)
            yadd = 0
            xOA = 0
            if srad < 0:
                xOA = -np.abs(srad)+np.sqrt(srad**2-mrad**2)
        elif surftype == 'parabolic':
            verts, faces, yadd, xOA, splitverts = sfc.add_parabolic_surface(-srad, mrad, N1, N2, theta, orig=opos,hole=self.cent_hole,hrad=hrad)
    nVerts = len(verts)
    
    #add side
    for j in range(N2):
        fi1 = nVerts+(j+1)%(N2)
        fi2 = nVerts+j
        fi3 = fi2-N2
        fi4 = fi1-N2
        faces.append([fi1,fi2,fi3,fi4])
    nVerts = len(verts)
    
    #add rear surface
    dvert, dfac, splitverts2 = sfc.add_flat_surface(mrad,N1,N2,-1,CT-xOA,yadd,nVerts=nVerts,hole=self.cent_hole,hrad=hrad)
    dvert = dvert[::-1]
    
    verts = verts+dvert
    faces = faces+dfac
        
    del dvert
    del dfac

    dummy1 = [0 for i in range(len(splitverts2))]
    dummy2 = [0 for i in range(len(splitverts))]
    splitverts = splitverts + dummy1
    splitverts2 = dummy2 + splitverts2[::-1]

    #fill hole
    if self.cent_hole:
        lv = len(verts)
        for j in range(N2):
            fi2 = j
            fi1 = (j+1)%N2
            fi4 = (j+1)%N2 + lv - N2
            fi3 = j + lv - N2
            faces.append([fi1,fi2,fi3,fi4])
    
    mesh = bpy.data.meshes.new(name="New Mirror")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)

    #custom split normals
    
    if not self.smooth_type:
        mesh = obj.data
        obj.select_set(True)
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [True, False, False]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        # mirror face
        for i in range(len(verts)):
            if splitverts[i] == 0:
                mesh.vertices[i].select=False
            else:
                mesh.vertices[i].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')

        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        #rear face
        for i in range(len(verts)):
            if splitverts2[i] == 0:
                mesh.vertices[i].select=False
            else:
                mesh.vertices[i].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')

        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        #hole
        if self.cent_hole:
            for i in range(N2):
                mesh.vertices[i].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')

            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

            nVerts = len(verts)
            for i in range(N2):
                mesh.vertices[nVerts -1 - i].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')

        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        #end split normal

    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()