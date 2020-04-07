import bpy
import numpy as np

from .. import surface as sfc
from .. import object_data_add
from .. import utils

def add_lens(self, context):
    edges = []
    
    srad1 = self.rad1
    srad2 = self.rad2
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
    
    #check surface radii for consistency
    ##check radius overflow
    if not utils.check_surface(np.abs(srad1), lrad): srad1=0
    if not utils.check_surface(np.abs(srad2), lrad): srad2=0
    ##check center thickness
    lsurf1, lsurf2 = 0, 0
    if not srad1 == 0:
        lsurf1 = srad1-ssig1*np.sqrt(srad1**2-lrad**2)
    if not srad2 == 0:
        lsurf2 = srad2-ssig2*np.sqrt(srad2**2-lrad**2)
    if (lsurf1 + lsurf2) > CT:
        CT = lsurf1 + lsurf2

    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces, splitverts = sfc.add_flat_surface(lrad,N1,N2, dshape=self.dshape)
    else:
        if self.ltype1 == 'spherical':
            verts, faces, splitverts = sfc.add_spherical_surface(srad1, lrad, N1, N2, dshape=self.dshape)
        elif self.ltype1 == 'aspheric':
            verts, faces, splitverts = sfc.add_aspheric_surface(srad1, k, A, lrad, N1, N2, dshape=self.dshape)
    
    
    nVerts = len(verts)
    
    #add surface2
    if srad2 == 0:
        #flat surface case
        dvert, dfac, splitverts2 = sfc.add_flat_surface(lrad,N1,N2,-1,CT,nVerts=nVerts, dshape=self.dshape)
        dvert = dvert[::-1]
    else:
        if self.ltype2 == 'spherical':
            dvert, dfac, splitverts2 = sfc.add_spherical_surface(srad2, lrad, N1, N2, -1, CT, nVerts=nVerts, dshape=self.dshape)
        elif self.ltype2 == 'aspheric':
            dvert, dfac, splitverts2 = sfc.add_aspheric_surface(srad2, k2, A2, lrad, N1, N2, -1, CT, nVerts=nVerts, dshape=self.dshape)
        #dvert, dfac = sfc.add_spherical_surface(srad2, lrad, N1, N2,-1, CT, nVerts=nVerts)
        dvert, dfac = dvert[::-1], dfac[::-1]
        
    verts = verts+dvert
    faces = faces+dfac

    nVerts2 = len(verts)

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

    dummy1 = [0 for i in range(len(splitverts2))]
    dummy2 = [0 for i in range(len(splitverts))]
    splitverts = splitverts + dummy1
    splitverts2 = dummy2 + splitverts2[::-1]
        
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

        for i in range(len(verts)):
            if splitverts2[i] == 0:
                mesh.vertices[i].select=False
            else:
                mesh.vertices[i].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')

        if self.dshape:
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

            for i in faces[-1]:
                mesh.vertices[i].select=True
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')

        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        #end split normals

    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()