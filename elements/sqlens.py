import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils

class OBJECT_OT_add_sqlens(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_sqlens"
    bl_label = "OptiCore Square Lens"
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
           default = 0,
           description="Radius of Curvature of Surface 2",
           unit = "LENGTH",
           )
    num1 : IntProperty(
           name="N",
           default = 32,
           description="Number of vertices (per axis)",
           min=3,
           )
    lenswidth : FloatProperty(
           name="Lens Width",
           default = 3.,
           description="Width of the lens",
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
    split_edge : BoolProperty(
            name="Edge Split",
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
        #layout.prop(self, 'ltype1')
        layout.prop(self, 'rad1')
        #layout.prop(self, 'ltype2')
        layout.prop(self, 'rad2')
        layout.prop(self, 'lenswidth')
        layout.prop(self, 'centerthickness')
        layout.prop(self, 'num1')
        #layout.prop(self, 'num2')
        if self.ltype1=='aspheric':
            layout.prop(self, 'k')
            layout.prop(self, 'A')
        if self.ltype2=='aspheric':
            layout.prop(self, 'k2')
            layout.prop(self, 'A2')
        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        layout.prop(self, 'shade_smooth')
        #layout.prop(self, 'split_edge')

    def execute(self, context):
        add_sqlens(self, context)
        return {'FINISHED'}

def add_sqlens(self, context):
    edges = []
    
    srad1 = self.rad1
    srad2 = self.rad2
    N1 = self.num1
    #N2 = self.num2
    N2 = N1
    lwidth = self.lenswidth
    CT = self.centerthickness
    #if self.ltype1 == 'aspheric':
    #    k = self.k
    #    A = self.A
    #if self.ltype2 == 'aspheric':
    #    k2 = self.k2
    #    A2 = self.A2

    ssig1 = 1
    if srad1 < 0:
        ssig1 = -1
    ssig2 = 1
    if srad2 < 0:
        ssig2 = -1
    
    #check surface radii for consistency
    ##check radius overflow
    if not utils.check_surface(np.abs(srad1), np.sqrt(2)/2*lwidth): srad1=0
    if not utils.check_surface(np.abs(srad2), np.sqrt(2)/2*lwidth): srad2=0
    ##check center thickness
    lsurf1, lsurf2 = 0, 0
    if not srad1 == 0:
        lsurf1 = srad1-ssig1*np.sqrt(srad1**2-0.5*lwidth**2)
    if not srad2 == 0:
        lsurf2 = srad2-ssig2*np.sqrt(srad2**2-0.5*lwidth**2)
    if (lsurf1 + lsurf2) > CT:
        CT = lsurf1 + lsurf2

    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces, splitverts = sfc.add_sqflat_surface(lwidth,N1,N2)
    else:
        if self.ltype1 == 'spherical':
            verts, faces, splitverts = sfc.add_sqspherical_surface(srad1, lwidth, N1, N2)
        #elif self.ltype1 == 'aspheric':
        #    verts, faces, splitverts = sfc.add_aspheric_surface(srad1, k, A, lrad, N1, N2)
    
    #add side
    nVerts = len(verts)

    
    if srad1==0 and srad2 == 0:
        faces.append([0,1] + [nVerts+2, nVerts+3])
        faces.append([2,3] + [nVerts, nVerts+1])
        faces.append([3,0] + [nVerts+3, nVerts])
        faces.append([1,2] + [nVerts+1, nVerts+2])
    elif srad1!=0 and srad2==0:
        faces.append([i for i in range(N2)] + [nVerts+2, nVerts+3])
        faces.append([i+N2*(N1-1) for i in range(N2)[::-1]] + [nVerts, nVerts+1])
        faces.append([i*N2 for i in range(N1)[::-1]] + [nVerts+3, nVerts])
        faces.append([i*N2 + N2-1 for i in range(N1)] + [nVerts+1, nVerts+2])
    elif srad1==0 and srad2!=0:
        faces.append([0,1] + [i+N2*(N1-1) + 4 for i in range(N2)])
        faces.append([2,3] + [i + 4 for i in range(N2)[::-1]])
        faces.append([3,0] + [i*N2 + N2-1 + 4 for i in range(N1)[::-1]])
        faces.append([1,2] + [i*N2 + 4 for i in range(N1)])
    else:
        faces.append([i for i in range(N2)] + [i+N2*(N1-1) + nVerts for i in range(N2)])
        faces.append([i+N2*(N1-1) for i in range(N2)[::-1]] + [i + nVerts for i in range(N2)[::-1]])
        faces.append([i*N2 for i in range(N1)[::-1]] + [i*N2 + N2-1 + nVerts for i in range(N1)[::-1]])
        faces.append([i*N2 + N2-1 for i in range(N1)] + [i*N2 + nVerts for i in range(N1)])
    
    #add surface2
    if srad2 == 0:
        #flat surface case
        dvert, dfac, splitverts2 = sfc.add_sqflat_surface(lwidth,N1,N2,-1,CT,nVerts=nVerts)
        dvert = dvert[::-1]
    else:
        if self.ltype2 == 'spherical':
            dvert, dfac, splitverts2 = sfc.add_sqspherical_surface(srad2, lwidth, N1, N2, -1, CT, nVerts=nVerts)
        #elif self.ltype2 == 'aspheric':
        #    dvert, dfac, splitverts2 = sfc.add_aspheric_surface(srad2, k2, A2, lrad, N1, N2, -1, CT, nVerts=nVerts)
        #dvert, dfac = sfc.add_spherical_surface(srad2, lrad, N1, N2,-1, CT, nVerts=nVerts)
        dvert, dfac = dvert[::-1], dfac[::-1]
        
    verts = verts+dvert
    faces = faces+dfac

    #dummy1 = [0 for i in range(len(splitverts2))]
    #dummy2 = [0 for i in range(len(splitverts))]
    #splitverts = splitverts + dummy1
    #splitverts2 = dummy2 + splitverts2[::-1]
    splitverts = splitverts + splitverts2
    
    del dvert
    del dfac

    mesh = bpy.data.meshes.new(name="New Square Lens")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)
    
    mesh = obj.data
    
    #custom split normals
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    bpy.ops.mesh.select_all(action='DESELECT')
    sel_mode = bpy.context.tool_settings.mesh_select_mode
    bpy.context.tool_settings.mesh_select_mode = [True, False, False]
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

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
    """
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
    """
    #end split normals

    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        bpy.ops.object.shade_smooth()
    if self.split_edge:
        obj.modifier_add(type='EDGE_SPLIT')