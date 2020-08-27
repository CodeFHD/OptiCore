import bpy
import numpy as np

from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from mathutils import Vector

from .. import surface as sfc
from .. import object_data_add
from .. import utils

class OBJECT_OT_add_CCretro(Operator, AddObjectHelper):
    """Cube Corner Retroreflector Array"""
    bl_idname = "mesh.add_ccretro"
    bl_label = "Retroreflector Array"
    bl_options = {'REGISTER', 'UNDO'}
    
    rshape : EnumProperty(
           name="Retroreflector Shape",
           items = {("circular","Circular",""),
                    ("rectangular","Rectangular","")},
           default = "rectangular",
           description="Outline Shape",
           #options={'HIDDEN'},
           )
    retrotype: EnumProperty(
           name="Retroreflector Type",
           items = {("cube","Cube Corner",""),
                    ("trirect","Trirectangular Tetrahedron","")},
           default = "cube",
           description="Reflector Type",
           #options={'HIDDEN'},
           )
    rrad : FloatProperty(
           name="Outline Radius",
           default = 12.,
           description="Radius of the Outline",
           unit = "LENGTH",
           )
    num1 : IntProperty(
           name="N",
           default = 64,
           description="Number of outline vertices",
           min=3,
           )
    rlenx : FloatProperty(
           name="Outline Length",
           default = 36.,
           description="Length of the Outline",
           unit = "LENGTH",
           )
    rleny : FloatProperty(
           name="Outline Width",
           default = 24.,
           description="Width of the Outline",
           unit = "LENGTH",
           )
    ccdist : FloatProperty(
           name="CubeCorner Spacing",
           default = 1.,
           description="Distance between CubeConer cnetres",
           unit = "LENGTH",
           )
    facethickness : FloatProperty(
           name="Face Thickness",
           default = -1.,
           description="thickness from CC to front/back face (Minimum 0.1% of CC depth)",
           unit = "LENGTH",
           )
    A : FloatVectorProperty(
           name="Tip Offset",
           default = (0.,0.,0.),
           description="XYZ offset of the cube tips",
           min=-1.,
           max=1.,
           )

    material_name : StringProperty(
            name="Material",
            default="",
           )

    shade_smooth : BoolProperty(
            name="Smooth Shading",
            default=False,
           )
    smooth_type : BoolProperty(
            name="Use Autosmooth",
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
        scene = context.scene
        #layout.prop(self, 'rshape')
        layout.prop(self, 'retrotype')
        if self.rshape == 'circular':
            layout.prop(self, 'rrad')
            layout.prop(self, 'num1')
        elif self.rshape == 'rectangular':
            layout.prop(self, 'rlenx')
            layout.prop(self, 'rleny')
        layout.prop(self, 'ccdist')
        layout.prop(self, 'facethickness')
        layout.prop(self, 'A')

        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')

    def execute(self, context):
        if self.retrotype == 'cube':
            add_ccretro_cube(self, context)
        elif self.retrotype == 'trirect':
            add_ccretro_trirect(self, context)
        return {'FINISHED'}



def add_ccretro_cube(self, context):
    """
    based on Trirectangular Tetrahedron
    """
    verts = []
    edges = []
    faces = []

    #define base parameters
    rrad = self.rrad
    num1 = self.num1
    rlenx = self.rlenx
    rleny = self.rleny
    ccdist = self.ccdist
    fthick = self.facethickness
    A = self.A

    #calculate base parameters
    h = ccdist/2
    dyrow = ccdist*np.sqrt(3)/2 #y offset betwee rows
    dxrow = ccdist/2 #x offset every second row
    rpts = ccdist/np.sqrt(3) #radius of corner points around center
    dz = ccdist/np.sqrt(6) #height offset of points, alternate dz and 2*dz

    ccdepth = 2*dz #total depth of cube corner

    A = [A[i]*ccdist/20 for i in range(len(A))]

    #calc number of cc elements:
    if self.rshape == 'circular':
        pass
    elif self.rshape == 'rectangular':
        nx = int(np.floor(rlenx/ccdist - 0.2))
        ny = int(np.floor(rleny/dyrow - 0.2))
        if ny%2 == 0:
            ny = ny-1
        ccxoff = -(nx-1)*ccdist/2
        ccyoff = -(ny-1)*dyrow/2

    #sanity checks
    if fthick < 0:
        facdir = -1
        if fthick > -1.01*ccdepth:
            fthick = -1.01*ccdepth
    else:
        facdir = 1
        if fthick < 0.01*ccdepth:
            fthick = 0.01*ccdepth

    #make outline verts
    if self.rshape == 'circular':
        pass
    elif self.rshape == 'rectangular':
        verts.append(Vector((-rlenx/2,-rleny/2,fthick)))
        verts.append(Vector((-rlenx/2,rleny/2,fthick)))
        verts.append(Vector((rlenx/2,rleny/2,fthick)))
        verts.append(Vector((rlenx/2,-rleny/2,fthick)))

        verts.append(Vector((-rlenx/2,-rleny/2,0)))
        verts.append(Vector((-rlenx/2,rleny/2,0)))
        verts.append(Vector((rlenx/2,rleny/2,0)))
        verts.append(Vector((rlenx/2,-rleny/2,0)))

        noutline = 8

    #make cc tip verts (in rows)
    for j in range(ny):
        for i in range(nx - j%2):
            x = ccxoff + i*ccdist + dxrow*(j%2) + A[0]
            y = ccyoff + j*dyrow + A[1]
            z = -ccdepth + A[2]
            verts.append(Vector((x,y,z)))
    ntips = ny*nx - ny//2

    #make cc top verts (in rows+1)
    for j in range(ny+1):
        for i in range(2*nx + 1):
            x = ccxoff + h*(i-1)
            y = ccyoff + j*dyrow - rpts/2*(1 + (i+j)%2)
            z = -1*((i+j-1)%2)*dz
            verts.append(Vector((x,y,z)))

    
    #fill faces
    
    #front and side faces
    faces.append([3,2,1,0])
    faces.append([0,1,5,4])
    faces.append([1,2,6,5])
    faces.append([2,3,7,6])
    faces.append([3,0,4,7])

    #interface
    #side1
    faces.append([7,4] + [noutline + ntips + 2*i + 1 for i in range(nx)])
    o = noutline + ntips
    for i in range(nx - 1):
        f1 = o + 2*i + 1
        f2 = o + 2*i + 2
        f3 = o + 2*i + 3
        faces.append([f1,f2,f3])
    #side2
    faces.append([4,5] + [noutline + ntips + (2*nx + 1)*(j+1) for j in range(ny)[::2]][::-1])
    for j in range(ny):
        if j%2 == 0:
            continue
        o = noutline + ntips + (2*nx + 1)*j
        f1 = o
        f2 = o + (2*nx + 1)*2
        f3 = o + (2*nx + 1) + 1
        faces.append([f1,f2,f3])
        
        f1 = o + 1
        f2 = o 
        f3 = o + (2*nx + 1) + 1
        faces.append([f1,f2,f3])
        
        f1 = o + (2*nx + 1)*2
        f2 = o + (2*nx + 1)
        f3 = o + (2*nx + 1) + 1
        faces.append([f1,f2,f3])
    #side3
    o = noutline + ntips + (2*nx + 1)*ny
    faces.append([5,6] + [o + 2*i for i in range(nx + 1)][::-1])
    for i in range(nx):
        o = noutline + ntips + (2*nx + 1)*ny + 2*i
        f1 = o
        f2 = o + 2
        f3 = o + 1
        faces.append([f1,f2,f3])
    #side4
    faces.append([6,7] + [noutline + ntips + (2*nx + 1)*(j+1) - 1 for j in range(ny + 1)[1::2]])
    for j in range(ny):
        if j%2 == 0:
            continue
        o = noutline + ntips + (2*nx + 1)*j
        f1 = o + (2*nx + 1)*2 -2
        f2 = o + (2*nx + 1)*3 -1
        f3 = o + (2*nx + 1) - 1
        faces.append([f1,f2,f3])
        
        f1 = o - 2
        f2 = o - 1
        f3 = o + (2*nx + 1) - 1
        faces.append([f1,f2,f3])
        
        f1 = o + (2*nx + 1)*2 - 2
        f2 = o + (2*nx + 1) - 1
        f3 = o + (2*nx + 1) - 2
        faces.append([f1,f2,f3])
    #close off corners
    #fl
    o = noutline + ntips
    f1 = o
    f2 = o + 1
    f3 = o + (2*nx + 1)
    faces.append([f1,f2,f3])
    f1 = 4
    f2 = o + (2*nx + 1)
    f3 = o + 1
    faces.append([f1,f2,f3])
    #fr
    o = noutline + ntips + (2*nx + 1)
    f1 = 7
    f2 = o - 2
    f3 = o + (2*nx + 1) - 1
    faces.append([f1,f2,f3])
    #rr
    o = noutline + ntips + (2*nx + 1)*(ny)
    f1 = o - 2
    f2 = o - 1
    f3 = o + (2*nx + 1) - 1
    faces.append([f1,f2,f3])

    #cubes
    for j in range(ny):
        for i in range(nx - j%2):
            o1 = noutline + j*nx - j//2
            f1 = o1 + i

            o2 = noutline + ntips + 2*i + j%2 + (2*nx + 1)*j
            o3 = o2 + 2*nx + 1
            
            f2 = o2 + 2
            f3 = o2 + 1
            f4 = o2
            faces.append([f1,f2,f3,f4])

            f2 = o3 + 1
            f3 = o3 + 2
            f4 = o2 + 2
            faces.append([f1,f2,f3,f4])
            
            f2 = o2
            f3 = o3
            f4 = o3 + 1
            faces.append([f1,f2,f3,f4])

            

    faces = [fac[::facdir] for fac in faces]

    

    mesh = bpy.data.meshes.new(name="New CC Retroreflector")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)
    
    """
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
    """

    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()












def add_ccretro_trirect(self, context):
    """
    based on Trirectangular Tetrahedron
    """
    verts = []
    edges = []
    faces = []

    #define base parameters
    rrad = self.rrad
    num1 = self.num1
    rlenx = self.rlenx
    rleny = self.rleny
    ccdist = self.ccdist
    fthick = self.facethickness
    A = self.A

    #calculate base parameters
    ccdepth = ccdist/np.sqrt(6) #depth of cube corner
    L = np.sqrt(3)/2*ccdist #spacing of the lines,; height of triangle as seen from top
    r = ccdist/(2*np.sqrt(3))
    yoffs = [r, L-r]

    A = [A[i]*ccdist/20 for i in range(len(A))]

    #calc number of cc elements:
    if self.rshape == 'circular':
        pass
    elif self.rshape == 'rectangular':
        nx = int(np.floor(rlenx/ccdist - 0.2))
        ny = int(np.floor(rleny/L - 0.2))
        if ny%2 == 1:
            ny = ny-1
        ccxoff = -nx*ccdist/2
        ccyoff = -ny*L/2

    #sanity checks
    if fthick < 0:
        facdir = -1
        if fthick > -1.01*ccdepth:
            fthick = -1.01*ccdepth
    else:
        facdir = 1
        if fthick < 0.01*ccdepth:
            fthick = 0.01*ccdepth

    #make outline verts
    if self.rshape == 'circular':
        pass
    elif self.rshape == 'rectangular':
        verts.append(Vector((rlenx/2,rleny/2,fthick)))
        verts.append(Vector((rlenx/2,-rleny/2,fthick)))
        verts.append(Vector((-rlenx/2,-rleny/2,fthick)))
        verts.append(Vector((-rlenx/2,rleny/2,fthick)))

        verts.append(Vector((rlenx/2,rleny/2,0)))
        verts.append(Vector((rlenx/2,-rleny/2,0)))
        verts.append(Vector((-rlenx/2,-rleny/2,0)))
        verts.append(Vector((-rlenx/2,rleny/2,0)))

        noutline = 8

    #make cc tip verts (in rows)
    for j in range(ny):
        for i in range(2*nx-1):
            x = ccxoff + (i+1)*ccdist/2 + A[0]
            y = ccyoff + j*L + yoffs[(j+i)%2] + A[1]
            z = -ccdepth + A[2]
            verts.append(Vector((x,y,z)))
    ntips = ny*(2*nx - 1)

    #make cc top verts (in rows+1)
    for j in range(ny+1):
        for i in range(nx + 1 - j%2):
            x = ccxoff + ccdist/2*(j%2) + ccdist*i
            y = ccyoff + j*L
            #z = 0
            verts.append(Vector((x,y,0)))


    #fill faces

    #front and side faces

    faces.append([3,2,1,0])
    faces.append([0,1,5,4])
    faces.append([1,2,6,5])
    faces.append([2,3,7,6])
    faces.append([3,0,4,7])

    #interface

    faces.append([5,6] + [noutline + ntips + i for i in range(nx+1)])
    faces.append([7,4] + [noutline + ntips + ny*(nx+1) - ny//2 + i for i in range(nx + 1 - ny%2)[::-1]])

    for j in range(int(ny//2)):
        o = noutline + ntips
        f1 = o + (2*nx + 1)*j
        f2 = o + (2*nx + 1)*j + 2*nx + 1
        f3 = o + (2*nx + 1)*j + nx + 1
        faces.append([f1,f2,f3])

        o = noutline + ntips + nx
        f1 = o + (2*nx + 1)*j
        f2 = o + (2*nx + 1)*j + nx
        f3 = o + (2*nx + 1)*j + 2*nx + 1
        faces.append([f1,f2,f3])

    o = noutline + ntips
    faces.append([6,7] + [o + i*(2*nx+1) for i in range(ny//2 + 1)[::-1]])
    o = noutline + ntips + nx
    faces.append([4,5] + [o + i*(2*nx+1) for i in range(ny//2 +1)])

    #cubes
    for j in range(ny):
        for i in range(2*nx-1):
            curo1 = (2*nx-1)*j + i
            f1 = noutline + curo1

            curo2 = (nx+1)*j - j//2

            if i%2 == 0 and j%2 == 0:
                f2 = noutline + ntips + curo2 + i//2 + 1
                f3 = noutline + ntips + curo2 + i//2
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + nx+1 + i//2
                f3 = noutline + ntips + curo2 + i//2 + 1
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + i//2
                f3 = noutline + ntips + curo2 + nx+1 + i//2
                faces.append([f1,f2,f3])

            if i%2 == 1 and j%2 == 0:
                f2 = noutline + ntips + curo2 + nx+1 + i//2 + 1
                f3 = noutline + ntips + curo2 + i//2 + 1
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + nx+1 + i//2
                f3 = noutline + ntips + curo2 + nx+1 + i//2 + 1
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + i//2 + 1
                f3 = noutline + ntips + curo2 + nx+1 + i//2
                faces.append([f1,f2,f3])



            if i%2 == 0 and j%2 == 1:
                f2 = noutline + ntips + curo2 + i//2
                f3 = noutline + ntips + curo2 + nx + i//2
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + nx + i//2
                f3 = noutline + ntips + curo2 + nx + i//2 + 1
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + nx + i//2 + 1
                f3 = noutline + ntips + curo2 + i//2
                faces.append([f1,f2,f3])

            if i%2 == 1 and j%2 == 1:
                f2 = noutline + ntips + curo2 + i//2 + 1
                f3 = noutline + ntips + curo2 + i//2
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + nx + i//2 + 1
                f3 = noutline + ntips + curo2 + i//2 + 1
                faces.append([f1,f2,f3])
                
                f2 = noutline + ntips + curo2 + i//2
                f3 = noutline + ntips + curo2 + nx + i//2 + 1
                faces.append([f1,f2,f3])

    faces = [fac[::facdir] for fac in faces]



    mesh = bpy.data.meshes.new(name="New CC Retroreflector")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)
    
    """
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
    """

    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        if self.smooth_type:
            obj.data.use_auto_smooth = 1
        bpy.ops.object.shade_smooth()