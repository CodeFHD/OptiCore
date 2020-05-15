import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import object_data_add
from .. import utils

from mathutils import Vector

import time

#brad = 0.01
#hrad = 0.003
#nph = 16
#nbev = 8



class OBJECT_OT_add_table(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_table"
    bl_label = "Breadboard"
    bl_options = {'REGISTER', 'UNDO'}
    
    length : FloatProperty(
           name="Length",
           default = 1.8,
           description="Length",
           unit = "LENGTH",
           )
    width : FloatProperty(
           name="Width",
           default = 0.9,
           description="Width",
           unit = "LENGTH",
           )
    thickness : FloatProperty(
           name="Thickness",
           default = 0.005,
           description="Thickness",
           unit = "LENGTH",
           )
    hrad : FloatProperty(
           name="Hole Radius",
           default = 0.003,
           description="Hole Radius",
           unit = "LENGTH",
           )
    nph : IntProperty(
           name="Verts per Hole",
           default = 16,
           description="Number of Vertices per hole (msut be divisible by 4)",
           min=4
           )
    hspac : FloatProperty(
           name="Hole Spacing",
           default = 0.025,
           description="Spacing between holes",
           unit = "LENGTH",
           )
    brad : FloatProperty(
           name="Bevel Radius",
           default = 0.003,
           description="Radius of edge bevel",
           unit = "LENGTH",
           )
    nbev : IntProperty(
           name="Bevel Verts",
           default = 6,
           description="Number of Vertices of the bevel.",
           min=3
           )
    shade_smooth : BoolProperty(
            name="Smooth Shading",
            default=True,
           )
    smooth_type : BoolProperty(
            name="Use Autosmooth (LuxCore v2.3)",
            default=True,
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

        #custom
        layout.prop(self, 'length')
        layout.prop(self, 'width')
        layout.prop(self, 'thickness')

        layout.prop(self, 'hrad')
        layout.prop(self, 'nph')
        layout.prop(self, 'hspac')
        layout.prop(self, 'brad')
        layout.prop(self, 'nbev')
        
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')

    def check_props(self):
        if self.length <= 0:
            self.length = 1.8
        if self.width <= 0:
            self.width = 0.9
        if self.thickness <= 0:
            self.length = 0.0127
        if self.hspac > 0.25*min(self.length, self.width):
            self.hspac = 0.25*min(self.length, self.width)
        if self.hrad >= self.hspac/2:
            self.hrad = 0.9*self.hspac/2
        if self.hrad <= 0:
            self.hrad = 0.01*self.hspac
        if self.brad <= 0:
            self.brad = 1e-5*min(self.length, self.width)


    def execute(self, context):
        self.check_props()
        add_table(self, context)
        return {'FINISHED'}

def add_table(self, context):

    leng = self.length
    wid = self.width
    thi = self.thickness
    hrad = self.hrad
    brad = self.brad
    hspac = self.hspac

    nph = self.nph
    if not nph%4 == 0:
        nph = nph//4
        nph = int(4*nph)
    nph2 = nph//2
    nph4 = nph//4
    nbev = self.nbev
    nbev4 = 4*nbev

    #general calc: number of holes per direction
    nx = int((leng-hspac)/hspac)
    ny = int((wid-hspac)/hspac)

    #first surface
    verts, edges, faces, splitverts = add_tableface(leng, wid, 0, nx, ny, hrad, brad, nph, nph4, nph2, nbev, nbev4, hspac, topface=True)
    nVerts = len(verts)
    nFaces = len(faces)
    
    #duplicate to bottom surface
    verts2, edges2, faces2, splitverts2 = add_tableface(leng, wid, -thi, nx, ny, hrad, brad, nph, nph4, nph2, nbev, nbev4, hspac)
    for i in range(len(faces2)):
        faces2[i] = [x + nVerts for x in faces2[i]]

    nVerts2 = len(verts2)
    nFaces2 = len(faces2)

    verts = verts + verts2
    edges = edges + edges2
    faces = faces + faces2

    #connect top and bottom
    ##sides
    for i in range(nbev4):
        faces.append([i,(i+1)%nbev4, (i+1)%nbev4 + nVerts, i + nVerts])
    ##holes
    o2 = nx*ny*nph
    for i in range(ny):
        for j in range(nx):
            o = nbev4 + (i*nx + j)*nph
            for k in range(nph):
                faces.append([o+o2+(k+1)%nph, o+o2+k, o+k+nVerts, o+(k+1)%nph + nVerts])

    nHoleFaces = nx*ny*nph
                
    #add object
    mesh = bpy.data.meshes.new(name="New Table")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    bpy.ops.object.mode_set(mode='OBJECT')

    #add split normals
    if not self.smooth_type:
        mesh = obj.data
        obj.select_set(True)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [True, False, False]
        bpy.ops.object.mode_set(mode='OBJECT')

        #edges top
        for i in range(nbev4):
            mesh.vertices[i].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        #edges bottom
        for i in range(nbev4):
            mesh.vertices[i+nVerts].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        #half of holes
        off1 = nbev4
        for i in range(ny):
            for j in range(nx):
                if i%2==0:
                    idx = 0
                else:
                    idx = 1
                idx2 = i*nx + j
                if idx2%2 == idx:
                    off2 = idx2*nph
                    for k in range(nph):
                        mesh.vertices[k+off1+off2].select=True
                else:
                    off2 = nVerts + idx2*nph
                    for k in range(nph):
                        mesh.vertices[k+off1+off2].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        #other half
        off1 = nbev4
        for i in range(ny):
            for j in range(nx):
                if i%2==0:
                    idx = 1
                else:
                    idx = 0
                idx2 = i*nx + j
                if idx2%2 == idx:
                    off2 = idx2*nph
                    for k in range(nph):
                        mesh.vertices[k+off1+off2].select=True
                else:
                    off2 = nVerts + idx2*nph
                    for k in range(nph):
                        mesh.vertices[k+off1+off2].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        bpy.ops.mesh.split_normals()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        #inside bevels
        o = nbev4 + nx*ny*nph
        for i in range(ny):
            print(i)
            for j in range(nx):
                idx = (i*nx + j)*nph
                for k in range(nph):
                    mesh.vertices[k + o + idx].select=True
            if i%8 == 0 or i == ny-1: #if all are computed at once blender crashes
                bpy.ops.object.mode_set(mode='EDIT', toggle=False)
                bpy.context.tool_settings.mesh_select_mode = sel_mode
                bpy.ops.mesh.split_normals()
                bpy.ops.mesh.select_all(action='DESELECT')
                bpy.ops.object.mode_set(mode='OBJECT')

    #select hole faces for easier assignment of hole material
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    sel_mode = bpy.context.tool_settings.mesh_select_mode
    bpy.context.tool_settings.mesh_select_mode = [False, False, True]
    bpy.ops.object.mode_set(mode='OBJECT')
    nfac = len(faces)
    for i in range(nx*ny*nph):
        mesh.polygons[nfac-i-1].select=True


    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()



def add_tableface(leng, wid, dz, nx, ny, hrad, brad, nph, nph4, nph2, nbev, nbev4, hspac, topface=False):
    verts, edges, faces, splitverts = [],[],[],[]

    #create outline verts
    ## four corners, spaced (L/2,W/2), 8 verts each
    xo = leng/2 - brad
    yo = wid/2 - brad
    xoffs = [-xo, -xo, xo, xo]
    yoffs = [-yo, yo, yo, -yo]
    for i in range(4):
        angoffset = (i+2)*np.pi/2
        xoff, yoff = xoffs[i], yoffs[i]
        for j in range(nbev):
            ang = j/(nbev-1)*np.pi/2 + angoffset
            verts.append(Vector((brad*np.sin(ang)+xoff,brad*np.cos(ang)+yoff,dz)))
            splitverts.append(1)

    #create hole-bevel verts
    if topface:
        for i in range(ny):
            yoff = (i/(ny-1)- 0.5)*(wid-hspac) 
            for j in range(nx):
                xoff = (j/(nx-1)- 0.5)*(leng-hspac)
                for k in range(nph):
                    ang = k/nph*2*np.pi + np.pi
                    verts.append(Vector((hrad*1.1*np.sin(ang)+xoff,hrad*1.1*np.cos(ang)+yoff,dz)))
                    splitverts.append(1)
    
    #create hole verts
    if topface:
        dz2 = -0.1*hrad
    else:
        dz2 = 0
    for i in range(ny):
        yoff = (i/(ny-1)- 0.5)*(wid-hspac) 
        for j in range(nx):
            xoff = (j/(nx-1)- 0.5)*(leng-hspac)
            for k in range(nph):
                ang = k/nph*2*np.pi + np.pi
                verts.append(Vector((hrad*np.sin(ang)+xoff,hrad*np.cos(ang)+yoff,dz+dz2)))
                splitverts.append(1)

    #surfaces four corners
    h1 = [0, nx*(ny-1)*nph + nph4, (nx*ny -1)*nph + 2*nph4, (nx-1)*nph + 3*nph4]
    h2 = [h1[0]+nph4, h1[1]+nph4, h1[2]+nph4, h1[3]-3*nph4]
    h1 = [x+4*nbev for x in h1]
    h2 = [x+4*nbev for x in h2]
    for i in range(4):
        fac = [x+i*nbev for x in range(nbev)] + [h2[i]] + [x+h1[i] for x in range(nph4)][::-1]
        faces.append(fac[::-1])

    #surfaces four edges
    for i in range(4):
        if i == 0:
            fac2 = [(i+1)%4*nph4 + j*nx*nph + 4*nbev for j in range(ny)][::-1]
        elif i == 1:
            fac2 = [(i+1)%4*nph4 + j*nph + nx*(ny-1)*nph + 4*nbev for j in range(nx)][::-1]
        elif i == 2:
            fac2 = [(i+1)%4*nph4 + (j+1)*nx*nph - nph + 4*nbev for j in range(ny)]
        else:
            fac2 = [(i+1)%4*nph4 + j*nph + 4*nbev for j in range(nx)]
        fac = [(i+1)*nbev - 1, ((i+1)*nbev)%(4*nbev)] + fac2
        faces.append(fac[::-1])
        
    #surfaces between holes
    dx = [0,1,1,0]
    dy = [0,0,1,1]
    for i in range(ny-1):
        for j in range(nx-1):
            fac = []
            for k in range(4):
                fac2 = [(nph2 - k*nph4 + x)%nph + ((i + dy[k])*nx + j + dx[k])*nph for x in range(nph4+1)]
                fac = fac + fac2
            fac = [x+4*nbev for x in fac]
            faces.append(fac)

    #surfaces hole-edge gaps
    ##left + right
    for i in range(ny-1):
        fac = [nph4 + x + i*nx*nph for x in range(nph4+1)] + [x + (i+1)*nx*nph for x in range(nph4+1)]
        fac = [x+4*nbev for x in fac]
        faces.append(fac)
        fac = [nph2 + x + ((i+1)*nx-1)*nph for x in range(nph4+1)] + [(3*nph4 + x)%nph + ((i+2)*nx-1)*nph for x in range(nph4+1)]
        fac = [x+4*nbev for x in fac]
        faces.append(fac)
    ##bottom and top
    for i in range(nx-1):
        fac = [(3*nph4 + x)%nph + i*nph for x in range(nph4+1)] + [x + i*nph + nph for x in range(nph4+1)]
        fac = [x+4*nbev for x in fac]
        faces.append(fac)
        fac = [nph2 + x + (nx*(ny-1) + i)*nph for x in range(nph4+1)] + [nph4 + x + (nx*(ny-1) + i)*nph + nph for x in range(nph4+1)]
        fac = [x+4*nbev for x in fac]
        faces.append(fac)

    #surfaces hole bevels
    if topface:
        nhVert = nx*ny*nph
        for i in range(ny):
            for j in range(nx):
                idx1 = (i*nx + j)*nph + nbev4
                idx2 = idx1 + nhVert
                for k in range(nph):
                    f1 = idx1 + k
                    f2 = idx1 + (k+1)%nph
                    f3 = idx2 + (k+1)%nph
                    f4 = idx2 + k
                    faces.append([f2,f1,f4,f3])

    return verts, edges, faces, splitverts