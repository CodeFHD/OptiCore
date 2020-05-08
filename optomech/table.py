import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import object_data_add
from .. import utils

from mathutils import Vector

brad = 0.01
hrad = 0.003
nph = 16
nbev = 8

if not nph%4 == 0:
    nph = nph//4
    nph = int(4*nph)
nph2 = nph//2
nph4 = nph//4
nbev4 = 4*nbev

class OBJECT_OT_add_table(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_table"
    bl_label = "OptiCore Optical Bench"
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
           default = 0.01,
           description="Thickness",
           unit = "LENGTH",
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

    def execute(self, context):
        add_table(self, context)
        return {'FINISHED'}

def add_table(self, context):

    leng = self.length
    wid = self.width
    thi = self.thickness

    #general calc: number of holes per direction
    nx = int((leng-0.025)/0.025)
    ny = int((wid-0.025)/0.025)

    #first surface
    verts, edges, faces = add_tableface(leng, wid, 0, nx, ny)
    nVerts = len(verts)
    
    #duplicate to bottom surface
    verts2, edges2, faces2 = add_tableface(leng, wid, -thi, nx, ny)
    for i in range(len(faces2)):
        faces2[i] = [x + nVerts for x in faces2[i]]

    verts = verts + verts2
    edges = edges + edges2
    faces = faces + faces2

    #connect top and bottom
    ##sides
    for i in range(nbev4):
        faces.append([i,(i+1)%nbev4, (i+1)%nbev4 + nVerts, i + nVerts])
    ##holes
    for i in range(ny):
        for j in range(nx):
            o = nbev4 + (i*nx + j)*nph
            for k in range(nph):
                faces.append([o+(k+1)%nph, o+k, o+k+nVerts, o+(k+1)%nph + nVerts])

    #add object
    mesh = bpy.data.meshes.new(name="New Table")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    mesh = obj.data
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='OBJECT')



def add_tableface(leng, wid, dz, nx, ny):
    verts, edges, faces = [],[],[]

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
    
    #create hole verts
    for i in range(ny):
        yoff = (i/(ny-1)- 0.5)*(wid-0.025) 
        for j in range(nx):
            xoff = (j/(nx-1)- 0.5)*(leng-0.025)
            for k in range(nph):
                ang = k/nph*2*np.pi + np.pi
                verts.append(Vector((hrad*np.sin(ang)+xoff,hrad*np.cos(ang)+yoff,dz)))

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

    return verts, edges, faces