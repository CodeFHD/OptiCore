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
           min=20,
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
            name="Use Autosmooth",
            default=True,
           )
    meshversion: BoolProperty(
            name="Mesh Type",
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
        #layout.prop(self, 'diam')
        layout.prop(self, 'length')
        layout.prop(self, 'nv')
        
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')
        #layout.prop(self, 'meshversion')

    def check_props(self):
        if self.length <= 20:
            self.length = 20


    def execute(self, context):
        self.check_props()
        if self.meshversion:
            add_post05_v2(self, context)
        else:
            add_post05_v1(self, context)
        return {'FINISHED'}

####################################################################################################
####################################################################################################
####################################################################################################

def add_post05_v1(self, context):
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
    njs = [int(n/2+1),n,int(n/2+1)]
    for j in range(3):
        verts2 = []
        nj = njs[j]
        for i in range(nj):
            ang = 2.*np.pi*i/n
            x = thrurad*np.cos(ang)
            y = np.sqrt(rad**2-x**2)
            if j == 0:
                z = 3.5
            elif j == 1:
                z = leng2+thrurad*np.sin(ang)
            else:
                z = leng-8
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

    ######################################
        
    #FACES
    oe2 = 2*int(n/2+1)
    oe = 4*int(n/2+1)
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
        fac = [i+4*n+oe2,(i+1)%n+4*n+oe2,(i+1)%n + 5*n+oe2,i+5*n+oe2]
        faces.append(fac[::-1])

    #top bevel with ring
    for j in range(6,11+1):
        for i in range(n):
            fac = [(i+1)%n+j*n+oe,i+j*n+oe,i+(j+1)*n+oe,(i+1)%n + (j+1)*n+oe]
            faces.append(fac)

    #top hole
    fac = [x+12*n+oe for x in range(n)]
    faces.append(fac[::-1])

    #holes to inner sides
    ##bottom
    for j in range(2):
        oj1 = 4*n if j%2 == 0 else 4*n + int(n/2+1)
        oj2 = 4*n + 2*int(n/2+1) if j%2 == 0 else 4*n + 2*int(n/2+1) + n
        fdir = 1 - 2*j
        for i in range(int(n/2)):
            fac = [i + oj1,
                   i + oj1 + 1,
                   (-i-1)%n + oj2,
                   (-i)%n + oj2]
            faces.append(fac[::fdir])
    ##top
    for j in range(2):
        oj1 = 0 if j%2 == 0 else n
        oj1 += 4*n + 2*int(n/2+1)
        oj2 = 0 if j%2 == 0 else int(n/2+1)
        oj2 += 4*n + 2*int(n/2+1)+ 2*n
        fdir = 1 - 2*j
        for i in range(int(n/2)):
            fac = [i + oj1,
                   i + oj1 + 1,
                   i + oj2 + 1,
                   i + oj2]
            faces.append(fac[::fdir])

    #sides of hole-to-side (weird name but ok)
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n)) + 1
    #
    f1 = 3*n + maxidx
    f2 = 4*n
    f3 = f2 + 2*int(n/2 + 1)
    f4 = f3 + 2*n
    f5 = f4 + 2*int(n/2 + 1) + maxidx
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac)
    #
    f1 = 3*n + n - maxidx
    f2 = 4*n + int(n/2)
    f3 = 4*n + 2*int(n/2 + 1) + int(n/2)
    f4 = f3 + 2*n
    f5 = f4 + int(n/2 + 1) + n - maxidx + 1
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac[::-1])
    #
    f1 = 3*n + int(n/2) - maxidx
    f2 = 4*n + int(n/2 + 1)
    f3 = 4*n + 2*int(n/2 + 1) + n
    f4 = 6*n + 2*int(n/2 + 1) + int(n/2) + 1
    f5 = 6*n + 4*int(n/2 + 1) + int(n/2) - maxidx
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac[::-1])
    #
    f1 = 3*n + int(n/2) + maxidx
    f2 = 4*n + 2*int(n/2 + 1) - 1
    f3 = 4*n + 2*int(n/2 + 1) + n + int(n/2)
    f4 = 6*n + 3*int(n/2 + 1) + int(n/2)
    f5 = 6*n + 4*int(n/2 + 1) + int(n/2) + maxidx
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac)
    
    #fill walls to holes
    #determine ang-range, determine outside wall vert-indices and close to side of hole
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n)) + 1
    ow1 = 3*n
    ow2 = 6*n + 4*int(n/2+1)
    oh11 = 4*n
    oh12 = 6*n + 2*int(n/2+1)
    oh21 = 4*n + int(n/2+1)
    oh22 = 6*n + 3*int(n/2+1)

    #per side, top and bottom half, loop through vert´s
    ## start vertex hole
    shs = [oh12,oh12,oh11,oh11] + [oh22,oh22,oh21,oh21]
    shs = [int(x) for x in shs]
    ##offset vertex hole
    ohs = [0,n/2,0,n/2] + [0,n/2,0,n/2]
    ohs = [int(x) for x in ohs]
    ## direction hole
    dirsh = [1,-1,1,-1] + [1,-1,1,-1]
    ## start vertex ring
    srs = [ow2,ow2,ow1,ow1] + [ow2,ow2,ow1,ow1]
    srs = [int(x) for x in srs]
    ## offset vertex ring
    ors = [0,0,0,0] + [n/2,n/2,n/2,n/2]
    ors = [int(x) for x in ors]
    ##direction ring
    dirsr = [-1,1,-1,1] + [1,-1,1,-1]
    #face directions
    facdirs = [1,-1,-1,1] + [-1,1,1,-1]
    #modulos
    #modulos = [n,n,-1,-1] + [n,n,-1,-1]
    modulos = [n,n,n,n] + [n,n,n,n]

    for j in range(8):
        shj = shs[j]
        ohj = ohs[j]
        dirh = dirsh[j]
        srj = srs[j]
        orj = ors[j]
        dirr = dirsr[j]
        fj = facdirs[j]
        modj = modulos[j]

        curidx = int(-1*maxidx*dirr)
        for i in range(int(n/4)):
            i *= dirh
            if modj > 0:
                dx1 = np.abs(verts[shj+(ohj+i)%modj][0] - verts[srj + (orj + curidx)%modj][0])
                dx2 = np.abs(verts[shj+(ohj+i)%modj][0] - verts[srj + (orj + curidx + dirr)%modj][0])
            else:
                dx1 = np.abs(verts[shj+(ohj+i)][0] - verts[srj + (orj + curidx)][0])
                dx2 = np.abs(verts[shj+(ohj+i)][0] - verts[srj + (orj + curidx + dirr)][0])
            if(dx2 < dx1):
                if modj > 0:
                    fac = [shj+(ohj+i)%modj, srj + (orj + curidx + dirr)%modj, srj + (orj + curidx)%modj]
                else:
                    fac = [shj+(ohj+i), srj + (orj + curidx + dirr), srj + (orj + curidx)]
                faces.append(fac[::fj])
                curidx += dirr
            if modj > 0:
                fac = [shj+(ohj+i)%modj, shj+(ohj+i+dirh)%modj, srj + (orj + curidx)%modj]
            else:
                fac = [shj+(ohj+i), shj+(ohj+i+dirh), srj + (orj + curidx)]
            faces.append(fac[::fj])
    

    #fill sides
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n)) + 1
    ow1 = 3*n
    ow2 = 6*n+4*int(n/2+1)
    for i in range(int(n/2-2*maxidx)):
        o = maxidx 
        fac = [ow1+o+i+1,ow1+o+i,ow2+o+i,ow2+o+i+1]
        faces.append(fac)
    for i in range(int(n/2-2*maxidx)):
        o = n/2 + maxidx
        fac = [ow1+o+i+1,ow1+o+i,ow2+o+i,ow2+o+i+1]
        faces.append(fac)

    mesh = bpy.data.meshes.new(name="New Post")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    bpy.ops.object.mode_set(mode='OBJECT')

    #custom split normals
    if not self.smooth_type:
        mesh = obj.data
        obj.select_set(True)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [True, False, False]
        bpy.ops.object.mode_set(mode='OBJECT')

        xo = 2*int(n/2+1)
        offsets = [0,
                   n,
                   2*n,
                   3*n,
                   4*n + xo,
                   5*n + xo,
                   6*n + 2*xo,
                   7*n + 2*xo,
                   8*n + 2*xo,
                   9*n + 2*xo,
                   10*n + 2*xo,
                   11*n + 2*xo,
                   12*n + 2*xo]

        for j in range(13):
            oj = offsets[j]
            for i in range(n):
                    mesh.vertices[oj+i].select=True
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT')

    if self.shade_smooth:
        if self.smooth_type:
            obj.data.auto_smooth_angle = 10*np.pi/180
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()

####################################################################################################
####################################################################################################
####################################################################################################

def add_post05_v2(self, context):
    verts, edges, faces, splitverts = [],[],[],[]

    rad = 12.7/2
    rad2 = rad-0.5
    rad3 = rad-1
    n = self.nv
    n2 = int(n/2)
    n4 = int(n/4)
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
    #compute angrange, create dummy ring of unity radius
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n))

    dummyverts = []

    for i in range(n4-maxidx):
        ang = i*2*np.pi/n
        dummyverts.append(Vector((np.cos(ang),np.sin(ang),0)))
    for i  in range(n2+1):
        ang = 2.*np.pi*i/n
        x = thrurad/rad*np.cos(ang)
        y = np.sqrt(1-x**2)
        ang = (i-n4)
        dummyverts.append(Vector((x,y,0)))
    for i in range(n2-2*maxidx-1):
        ang = (i-n4+maxidx+1)*2*np.pi/n + np.pi
        dummyverts.append(Vector((np.cos(ang),np.sin(ang),0)))
    for i  in range(n2+1):
        ang = 2.*np.pi*i/n
        x = -1.*thrurad/rad*np.cos(ang)
        y = -1*np.sqrt(1-x**2)
        ang = (i-n4)
        dummyverts.append(Vector((x,y,0)))
    for i in range(n4-maxidx - 1)[::-1]:
        ang = -(i+1)*2*np.pi/n
        dummyverts.append(Vector((np.cos(ang),np.sin(ang),0)))

    npr = len(dummyverts) #n_per_ring

    #add thru holes first
    verts2 = []
    for i in range(n):
        ang = 2.*np.pi*i/n
        x = thrurad*np.cos(ang)
        y = np.sqrt(rad**2-x**2)
        z = leng2+thrurad*np.sin(ang)
        verts.append(Vector((x,y,z)))
        verts2.append(Vector((x,-y,z)))
    for v in verts2:
        verts.append(v)
    del verts2

    #add other rings
    rrads = [3.0,3.0,
             rad2,rad,rad,
             rad - 1*2/5, rad - 1*2.5/5 - .5, rad - 1*3/5,
             rad3, 2.0, 2.0]
    rzs = [8.6, 0,
           0, 0.5, leng-5,
           leng - 3, leng - 2.5, leng - 2,
           leng, leng, leng - 7.]

    for crad, z in np.column_stack((rrads, rzs)):
        for v in dummyverts:
            verts.append(Vector((crad*v[0],crad*v[1],z)))

    #FACES

    #thru hole
    for i in range(n):
        fac = [(i+1)%n, i, i+n, (i+1)%n + n]
        faces.append(fac)

    #bottom of M6 hole
    fac = [x + 2*n for x in range(npr)]
    faces.append(fac[::-1])

    #bottom part
    for j in range(3):
        o = 2*n + j*npr
        for i in range(npr):
            fac = [i + o, (i+1)%npr + o, (i+1)%npr + npr + o, i + npr + o]
            faces.append(fac)

    #sides to holes
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n))
    
    ##h1 bot
    o = 2*n + 3*npr  + n4 - maxidx#(n4 - maxidx) + (n2+1) + 2*(n4-maxidx)
    for k in range(n2):
        i = k + o
        j = n - k - 1
        fac = [i, i+1, j, (j+1)%n]
        faces.append(fac)
    ##h1 top
    o = 2*n + 4*npr + n4 - maxidx
    for k in range(n2):
        i = k + o
        j = k
        fac = [i, i+1, (j+1)%n, j]
        faces.append(fac[::-1])
    ##h2 bot
    o = 2*n + 3*npr  + (n4 - maxidx) + 2*(n2+1) + 2*(n4-maxidx) - 3
    for k in range(n2):
        i = o - k
        j = n - k - 1
        fac = [i, i+1, (j+1)%n + n, j + n]
        faces.append(fac)
    ##h2 top
    o =  2*n + 4*npr  + (n4 - maxidx) + 2*(n2+1) + 2*(n4-maxidx) - 3
    for k in range(n2):
        i = o - k
        j = k
        fac = [i, i+1, j + n, (j+1)%n + n]
        faces.append(fac[::-1])

    ##non-hole sides
    angrange = np.arcsin(thrurad/rad)
    maxidx = int(angrange//(2*np.pi/n)) + 1

    o = 2*n + 3*npr
    for j in range(n2-2*maxidx):
        i = (j - (n4 - maxidx))%npr
        fac = [i + o, (i+1)%npr + o, (i+1)%npr + o + npr, i + o + npr]
        faces.append(fac)
    for j in range(n2-2*maxidx):
        i = (j - (n4 - maxidx) + int(npr/2))%npr
        fac = [i + o, (i+1)%npr + o, (i+1)%npr + o + npr, i + o + npr]
        faces.append(fac)

    ##interconnect
    ##side1
    o = 2*n + 3*npr

    f1 = 0
    f2 = o + n4 - maxidx + 1
    f3 = f2 - 1
    f4 = o + npr + n4 - maxidx
    f5 = f4 + 1
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac[::-1])
    #
    f1 = n2
    f2 = o + n4 - maxidx + 1 + n2
    f3 = f2 + 1 
    f4 = o + npr + n4 - maxidx  + n2 + 2
    f5 = f4 - 1
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac)

    ##side2
    o = 2*n + 3*npr + 2*(n4-maxidx + 1) + n2

    f1 = n2 + n
    f2 = o + n4 - maxidx + 1
    f3 = f2 - 1
    f4 = o + npr + n4 - maxidx
    f5 = f4 + 1
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac[::-1])
    #
    f1 = 0 + n
    f2 = o + n4 - maxidx + 1 + n2
    f3 = f2 + 1 
    f4 = o + npr + n4 - maxidx  + n2 + 2
    f5 = f4 - 1
    fac = [f1,f2,f3,f4,f5]
    faces.append(fac)

    #top part
    for j in range(6):
        o = 2*n + (j+4)*npr
        for i in range(npr):
            fac = [i + o, (i+1)%npr + o, (i+1)%npr + npr + o, i + npr + o]
            faces.append(fac)

    #bottom of M4 hole
    fac = [x + 2*n + 10*npr for x in range(npr)]
    faces.append(fac)

    mesh = bpy.data.meshes.new(name="New Post")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    bpy.ops.object.mode_set(mode='OBJECT')

    #custom split normals
    if not self.smooth_type:
        mesh = obj.data
        obj.select_set(True)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [True, False, False]
        bpy.ops.object.mode_set(mode='OBJECT')

        offsets = [0,
                   1*n,
                   2*n,
                   2*n + 1*npr,
                   2*n + 2*npr,
                   2*n + 3*npr,
                   2*n + 4*npr,
                   2*n + 5*npr,
                   2*n + 6*npr,
                   2*n + 7*npr,
                   2*n + 8*npr,
                   2*n + 9*npr,
                   2*n + 10*npr]

        for j in range(13):
            oj = offsets[j]
            if j < 2:
                nrange = n
            else:
                nrange = npr
            for i in range(nrange):
                    mesh.vertices[oj+i].select=True
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.context.tool_settings.mesh_select_mode = sel_mode
            bpy.ops.mesh.split_normals()
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT')

    if self.shade_smooth:
        if self.smooth_type:
            obj.data.auto_smooth_angle = 10*np.pi/180
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()