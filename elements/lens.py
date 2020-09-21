import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils

class OBJECT_OT_add_lens(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "Lens"
    bl_options = {'REGISTER', 'UNDO'}
    
    makedoublet : BoolProperty(
            name="Doublet Lens",
            default=False,
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
    rad3 : FloatProperty(
           name="Surface 3 Radius",
           default = 24.,
           description="Radius of Curvature of Surface 3",
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
    centerthickness2 : FloatProperty(
           name="Center Thickness 2",
           default = 1.,
           description="Center thickness of lens segment 2",
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
    k3 : FloatProperty(
           name="k3",
           default = 0.,
           description="Aspheric conical constant",
           )
    A3 : FloatVectorProperty(
           name="A3",
           default = (0.,0.,0.),
           description="Aspheric correction coefficients",
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
            name="Use Autosmooth",
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
        if self.makedoublet:
            layout.prop(self, 'ltype3')
            layout.prop(self, 'rad3')
        layout.prop(self, 'lensradius')
        layout.prop(self, 'centerthickness')
        if self.makedoublet:
            layout.prop(self, 'centerthickness2')
        layout.prop(self, 'num1')
        layout.prop(self, 'num2')
        if self.ltype1=='aspheric':
            layout.prop(self, 'k')
            layout.prop(self, 'A')
        if self.ltype2=='aspheric':
            layout.prop(self, 'k2')
            layout.prop(self, 'A2')
        if self.makedoublet:
            if self.ltype3=='aspheric':
                layout.prop(self, 'k3')
                layout.prop(self, 'A3')
        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        if self.makedoublet:
            layout.prop_search(self, "material_name2", bpy.data, "materials", icon="NONE")
            layout.prop_search(self, "material_name3", bpy.data, "materials", icon="NONE")
        layout.prop(self, 'makedoublet')
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')
        layout.prop(self, 'dshape')
        layout.prop(self, 'optiverts')

    def execute(self, context):
        add_lens(self, context)
        return {'FINISHED'}

def add_lens(self, context):
    edges = []
    
    srad1 = self.rad1
    srad2 = -self.rad2
    N1 = self.num1
    N2 = self.num2
    lrad = self.lensradius
    CT = self.centerthickness
    if self.ltype1 == 'aspheric':
        k = self.k
        A = self.A
    if self.ltype2 == 'aspheric':
        k2 = self.k2
        A2 = self.A2
    ssig1 = 1
    if srad1 < 0:
        ssig1 = -1
    ssig2 = 1
    if srad2 < 0:
        ssig2 = -1

    md = self.makedoublet
    if md:
        srad3 = -self.rad3
        CT2 = self.centerthickness2
        if self.ltype3 == 'aspheric':
            k3 = self.k3
            A3 = self.A3
        ssig3 = 1
        if srad3 < 0:
            ssig3 = -1
    
    #check surface radii for consistency
    ##check radius overflow
    if not utils.check_surface(np.abs(srad1), lrad): srad1=0
    if not utils.check_surface(np.abs(srad2), lrad): srad2=0
    if md:
        if not utils.check_surface(np.abs(srad3), lrad): srad3=0
    ##check center thickness
    lsurf1, lsurf2 = 0, 0
    if not srad1 == 0:
        lsurf1 = srad1-ssig1*np.sqrt(srad1**2-lrad**2)
    if not srad2 == 0:
        lsurf2 = srad2-ssig2*np.sqrt(srad2**2-lrad**2)
    if (lsurf1 + lsurf2) > CT:
        CT = lsurf1 + lsurf2
    if md:
        lsurf3 = 0
        if not srad3 == 0:
            lsurf3 = srad3-ssig3*np.sqrt(srad3**2-lrad**2)
        if (-lsurf2 + lsurf3) > CT2:
            CT2 = -lsurf2 + lsurf3

    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces, splitverts = sfc.add_flat_surface(lrad,N1,N2, dshape=self.dshape)
    else:
        if self.ltype1 == 'spherical':
            verts, faces, splitverts, N1, N2 = sfc.add_spherical_surface(srad1, lrad, N1, N2, dshape=self.dshape, optiverts=self.optiverts)
        elif self.ltype1 == 'aspheric':
            verts, faces, splitverts = sfc.add_aspheric_surface(srad1, k, A, lrad, N1, N2, dshape=self.dshape)
    
    
    nVerts = len(verts)
    nFacs = len(faces)
    
    #add surface2
    if srad2 == 0:
        #flat surface case
        dvert, dfac, splitverts2 = sfc.add_flat_surface(lrad,N1,N2,-1,CT,nVerts=nVerts, dshape=self.dshape)
        dvert = dvert[::-1]
    else:
        if self.ltype2 == 'spherical':
            dvert, dfac, splitverts2, N12, N22 = sfc.add_spherical_surface(srad2, lrad, N1, N2, -1, CT, nVerts=nVerts, dshape=self.dshape, optiverts=self.optiverts)
        elif self.ltype2 == 'aspheric':
            dvert, dfac, splitverts2 = sfc.add_aspheric_surface(srad2, k2, A2, lrad, N1, N2, -1, CT, nVerts=nVerts, dshape=self.dshape)
        #dvert, dfac = sfc.add_spherical_surface(srad2, lrad, N1, N2,-1, CT, nVerts=nVerts)
        dvert, dfac = dvert[::-1], dfac[::-1]

    #reverse face normals if doublet lens
    if md:
        dfac = [f[::-1] for f in dfac]
        
    verts = verts+dvert
    faces = faces+dfac

    nVerts2 = len(verts)
    nFacs2 = len(faces)
    ndFac2 = len(dfac)

    #add side
    for j in range(N2-self.dshape):
        fi1 = nVerts+(j+1)%(N2)
        fi2 = nVerts+j
        fi3 = fi2-N2
        fi4 = fi1-N2
        faces.append([fi1,fi2,fi3,fi4])
    if self.dshape:
        if srad1==0:
            ptss1 = [N2-1,0]
        else:
            ptss1 = [(x+1)*N2 for x in range(N1)[::-1]] + [0] + [x*N2+1 for x in range(N1)]
        if srad2==0:
            ptss2 = [nVerts2-N2, nVerts2-1]
        else:
            ptss2 = [nVerts + x*N2 for x in range(N1)] + [nVerts2-1] + [nVerts + (x+1)*N2-1 for x in range(N1)[::-1]]

        faces.append(ptss1 + ptss2)

    nFacSide = N2

    ####################################
    if md:
        #add surface3
        if srad3 == 0:
            #flat surface case
            dvert, dfac, splitverts3 = sfc.add_flat_surface(lrad,N1,N2,-1,CT+CT2,nVerts=nVerts2, dshape=self.dshape)
            dvert = dvert[::-1]
        else:
            if self.ltype3 == 'spherical':
                dvert, dfac, splitverts3, N12, N22 = sfc.add_spherical_surface(srad3, lrad, N1, N2, -1, CT+CT2, nVerts=nVerts2, dshape=self.dshape, optiverts=self.optiverts)
            elif self.ltype3 == 'aspheric':
                dvert, dfac, splitverts3 = sfc.add_aspheric_surface(srad3, k2, A2, lrad, N1, N2, -1, CT+CT2, nVerts=nVerts2, dshape=self.dshape)
            #dvert, dfac = sfc.add_spherical_surface(srad2, lrad, N1, N2,-1, CT, nVerts=nVerts)
            dvert, dfac = dvert[::-1], dfac[::-1]
        
        verts = verts+dvert
        faces = faces+dfac

        nvs1 = N2*N1+1
        if srad1 == 0:
            nvs1 = N2
        nvs2 = N2*N1+1
        if srad2 == 0:
            nvs2 = N2

        nVerts3 = len(verts)
        nFacs3 = len(faces)
        ndFac3 = len(dfac)

        #add side
        for j in range(N2-self.dshape):
            fi1 = nVerts2+(j+1)%(N2)
            fi2 = nVerts2+j
            fi3 = fi2-nvs2
            fi4 = fi1-nvs2
            faces.append([fi1,fi2,fi3,fi4])
        
        if self.dshape:
            if srad2==0:
                ptss1 = [nvs1 + N2-1, nvs1]
            else:
                ptss1 = [nVerts + (x+1)*N2-1 for x in range(N1)] + [nVerts2-1] + [nVerts + x*N2 for x in range(N1)[::-1]]   #[nvs1 + (x+1)*N2 for x in range(N1)[::-1]] + [nvs1] + [nvs1 + x*N2+1-N2 for x in range(N1)]
            if srad3==0:
                ptss2 = [nVerts3-N2, nVerts3-1]
            else:
                ptss2 = [nVerts2 + x*N2 for x in range(N1)] + [nVerts3-1] + [nVerts2 + (x+1)*N2-1 for x in range(N1)[::-1]]
        
            faces.append(ptss1 + ptss2)

    ####################################

    dummy1 = [0 for i in range(len(splitverts))]
    dummy2 = [0 for i in range(len(splitverts2))]
    if md:
        dummy3 = [0 for i in range(len(splitverts3))]
    else:
        dummy3 = []
    splitverts = splitverts + dummy2 + dummy3
    splitverts2 = dummy1 + splitverts2[::-1] + dummy3
    if md:
        splitverts3 = dummy1 + dummy2 + splitverts3[::-1]

    #dummy1 = [0 for i in range(len(splitverts2))]
    #dummy2 = [0 for i in range(len(splitverts))]
    #splitverts = splitverts + dummy1
    #splitverts2 = dummy2 + splitverts2[::-1]
        
    del dvert
    del dfac
    mesh = bpy.data.meshes.new(name="New Lens")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)
    if not self.smooth_type:
        mesh = obj.data
    
        #custom split normals
        obj.select_set(True)
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [True, False, False]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        #splitverts 1
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

        #splitverts 2
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

        #splitverts 
        if md:
            for i in range(len(verts)):
                if splitverts3[i] == 0:
                    mesh.vertices[i].select=False
                else:
                    mesh.vertices[i].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        #dshape splits
        if self.dshape:
            for i in faces[-1]:
                mesh.vertices[i].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            if md:
                for i in faces[-ndFac3 - nFacSide - 1]:
                    mesh.vertices[i].select=True
                bpy.ops.object.mode_set(mode='EDIT', toggle=False)
                bpy.context.tool_settings.mesh_select_mode = sel_mode
                bpy.ops.mesh.split_normals()
                bpy.ops.mesh.select_all(action='DESELECT')
                bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        
        #end split normals

    #assign material(s)
    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)

    if md:
        cond2 = self.material_name2 in bpy.data.materials
        cond3 = self.material_name3 in bpy.data.materials
    if md and cond2 and cond3:
        if self.material_name2 in bpy.data.materials:
            mat2 = bpy.data.materials[self.material_name2]
            obj.data.materials.append(mat2)
        if self.material_name3 in bpy.data.materials:
            mat3 = bpy.data.materials[self.material_name3]
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

    #apply smooth shading
    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()