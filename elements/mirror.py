import bpy
import numpy as np

from .. import surface as sfc
from .. import object_data_add
from .. import utils

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
    mesh = obj.data

    #custom split normals
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
    #end split normals


    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        bpy.ops.object.shade_smooth()
    if self.split_edge:
        bpy.ops.object.modifier_add(type='EDGE_SPLIT')