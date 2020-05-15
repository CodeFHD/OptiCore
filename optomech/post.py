import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import object_data_add
from .. import utils

from mathutils import Vector

import time

class OBJECT_OT_add_post(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_post"
    bl_label = "Post"
    bl_options = {'REGISTER', 'UNDO'}
    
    diam : EnumProperty(
           name="Diameter",
           items = {("0.5","0.5-inch",""),
                    ("1.0","1.0-inch",""),
                    ("1.5","1.0-inch","")},
           default = "0.5",
           description="Post diameter",
           #options={'HIDDEN'},
           )
    length : FloatProperty(
           name="Length",
           default = 30,
           description="Length",
           unit = "LENGTH",
           )
    nv : IntProperty(
           name="Number of Vertices",
           default = 64,
           description="Number of Vertices",
           min=16
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
        layout.prop(self, 'diam')
        layout.prop(self, 'length')
        layout.prop(self, 'nv')
        
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')

    def check_props(self):
        if self.length <= 5:
            self.length = 5


    def execute(self, context):
        self.check_props()
        add_post05(self, context)
        return {'FINISHED'}
    
def add_post05(self, context):
    verts, edges, faces, splitverts = [],[],[],[]

    rad = 12.7/2
    rad2 = rad-0.5
    rad3 = rad-1
    n = self.nv
    leng = self.length
    if leng < 25:
        dl = 9.5
    elif leng < 150:
        dl = 10.2
    else:
        dl = 12.7
    leng2 = leng-dl
    thrurad = 1.6

    #VERTS
    #bottom hole M6
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((3.0*np.sin(ang),3.0*np.cos(ang),8.6)))
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((3.0*np.sin(ang),3.0*np.cos(ang),0)))

    #base with bevel
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((rad2*np.sin(ang),rad2*np.cos(ang),0)))

    #base at full diameter
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((rad*np.sin(ang),rad*np.cos(ang),0.5)))

    #Thru hole on side
    verts2 = []
    for i in range(n):
        ang = 2.*np.pi*i/n
        x = thrurad*np.cos(ang)
        z = leng2+thrurad*np.sin(ang)
        y = np.sqrt(rad**2-x**2)
        verts.append(Vector((x,y,z)))
        verts2.append(Vector((x,-y,z)))
    for v in verts2:
        verts.append(v)
    del verts2


    #end of full diameter
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((rad*np.sin(ang),rad*np.cos(ang),leng-5)))


    #indentring
    radx = rad - 1*2/5
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((radx*np.sin(ang),radx*np.cos(ang),leng-3)))
    radx = rad - 1*2.5/5 - .5
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((radx*np.sin(ang),radx*np.cos(ang),leng-2.5)))
    radx = rad - 1*3/5
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((radx*np.sin(ang),radx*np.cos(ang),leng-2)))

    #top
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((rad3*np.sin(ang),rad3*np.cos(ang),leng)))


    #top hole M4
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((2.*np.sin(ang),2.*np.cos(ang),leng)))
    for i in range(n):
        ang = 2.*np.pi*i/n
        verts.append(Vector((2.*np.sin(ang),2.*np.cos(ang),leng-7.)))

    #FACES
    #bottom hole
    fac = [x for x in range(n)]
    faces.append(fac)
    for i in range(n):
        fac = [(i+1)%n,i,i+n,(i+1)%n + n]
        faces.append(fac)

    #bottom face
    for i in range(n):
        fac = [(i+1)%n+n,i+n,i+2*n,(i+1)%n + 2*n]
        faces.append(fac)

    #bottom bevel
    for i in range(n):
        fac = [(i+1)%n+2*n,i+2*n,i+3*n,(i+1)%n + 3*n]
        faces.append(fac)

    #Thru hole
    for i in range(n):
        fac = [i+4*n,(i+1)%n+4*n,(i+1)%n + 5*n,i+5*n]
        faces.append(fac)

    #top bevel with ring
    for j in range(6,11+1):
        for i in range(n):
            fac = [(i+1)%n+j*n,i+j*n,i+(j+1)*n,(i+1)%n + (j+1)*n]
            faces.append(fac)

    #top hole
    fac = [x+12*n for x in range(n)]
    faces.append(fac[::-1])
    #fill walls to holes
    #determine ang-range, determine outside wall vert-indices and close to side of hole
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n)) + 1
    ow1 = 3*n
    oh1 = 4*n
    oh2 = 5*n
    ow2 = 6*n

    fac = [ow1 + maxidx, oh1, ow2 + maxidx]
    faces.append(fac)
    fac = [ow1 + n - maxidx, oh1 + n/2, ow2 + n - maxidx]
    faces.append(fac)
    fac = [ow1 + n/2 + maxidx, oh2 + n/2, ow2 + n/2 + maxidx]
    faces.append(fac)
    fac = [ow1 + n/2 - maxidx, oh2, ow2 + n/2 - maxidx]
    faces.append(fac)
    #per side, top and bottom half, loop through vert´s
    #determine sh*t and do it, i.e. TODO documentation
    ohcs = [oh1,oh1+n/2,oh1,oh1+n/2,oh2+n/2,oh2,oh2+3*n/4,oh2+n/2]
    ohcs = [int(x) for x in ohcs]
    #owcs = [ow1,ow1,ow2,ow2,ow1+n/2,ow1+n/2,ow2+n/2,ow2+n/2]
    owcs = [ow2,ow2,ow1,ow1,ow2+n/2,ow2+n/2,ow1+n/2,ow1+n/2]
    owcs = [int(x) for x in owcs]
    dirs = [1,-1,-1,1,-1,1,-1,1]
    dirs2 = [1,-1,1,-1,1,-1,1,-1]
    dirs3 = [0,n/2,0,n/2,0,n/2,0,n/2]
    dirs3 = [int(x) for x in dirs3]

    for j in range(8):
        ohc = ohcs[j]
        owc = owcs[j]
        dic = dirs[j]
        dic2 = dirs2[j]
        dic3 = dirs3[j]
        curidx = dic2*maxidx
        for i in range(int(n/4)):
            i += dic3
            dx1 = np.abs(verts[ohc+i%n][0] - verts[owc+curidx%n][0])
            dx2 = np.abs(verts[ohc+i%n][0] - verts[owc+(curidx-dic)%n][0])
            if dic*dic2 < 0:
                dx1,dx2 = dx2,dx1
            if(dx2 < dx1):
                fac = [ohc+i%n, owc + (curidx-dic)%n, owc + curidx%n]
                faces.append(fac[::dic])
                curidx -= dic2
            fac = [ohc+i%n, ohc + (i+dic)%n, owc + curidx%n]
            faces.append(fac[::dic])


    #custom split normals
    ##TODO loop through all rings of n verts and split

    mesh = bpy.data.meshes.new(name="New Post")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    bpy.ops.object.mode_set(mode='OBJECT')

    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()