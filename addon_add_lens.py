bl_info = {
    "name": "LensCore",
    "author": "Johannes Hinrichs (CodeHD)",
    "version": (1, 0),
    "blender": (2, 79, 0),
    "location": "View3D > Add > Mesh > Add Lens",
    "description": "Adds a new optical lens or mirror",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh",
    }


import bpy
import numpy as np
from bpy.types import Operator
from bpy.props import FloatVectorProperty, FloatProperty, IntProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector

def add_lens(self, context):
    edges = []
    
    srad1 = self.rad1
    srad2 = self.rad2
    N1 = self.num1
    N2 = self.num2
    lrad = self.lensradius
    CT = self.centerthickness
    
    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces = add_flat_surface(lrad,N1,N2)
    else:
        verts, faces = add_curved_surface(srad1, lrad, N1, N2)
    nVerts = len(verts)
    
    #add side
    for j in range(N2):
        fi1 = nVerts+(j+1)%(N2)
        fi2 = nVerts+j
        fi3 = fi2-N2
        fi4 = fi1-N2
        faces.append([fi1,fi2,fi3,fi4])
    nVerts = len(verts)
    
    #add surface2
    if srad2 == 0:
        #flat surface case
        dvert, dfac = add_flat_surface(lrad,N1,N2,-1,CT,nVerts)
        dvert = dvert[::-1]
    else:
        dvert, dfac = add_curved_surface(srad2, lrad, N1, N2,-1, CT, nVerts)
        dvert, dfac = dvert[::-1], dfac[::-1]
        
    verts = verts+dvert
    faces = faces+dfac
        
    del dvert
    del dfac

    mesh = bpy.data.meshes.new(name="New Lens")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    object_data_add(context, mesh, operator=self)

def add_mirror(self, context):
    edges = []
    
    srad = self.rad
    N1 = self.num1
    N2 = self.num2
    mrad = self.mirrorradius
    CT = self.centerthickness
    
    #compute mirror surface
    verts, faces = add_curved_surface(-1.*srad, mrad, N1, N2)
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
    dvert, dfac = add_flat_surface(mrad,N1,N2,-1,CT,nVerts)
    dvert = dvert[::-1]
    
    verts = verts+dvert
    faces = faces+dfac
        
    del dvert
    del dfac
    
    mesh = bpy.data.meshes.new(name="New Mirror")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    object_data_add(context, mesh, operator=self)


def add_flat_surface(lrad,N1,N2,nsurf=1,xadd=0,nVerts=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """
    
    surfadd=0
    if nsurf == -1:
        surfadd = nVerts+N2-1
    
    verts = []
    faces = []
    
    for j in range(N2)[::nsurf]:
        b = 2*np.pi*j/N2
        verts.append(Vector((-1.*xadd,lrad*np.sin(b),lrad*np.cos(b))))
    faces.append([int(surfadd+nsurf*x) for x in range(N2)[::-1]])
    
    return verts, faces

def add_curved_surface(rad,lrad,N1,N2,nsurf=1,xadd=0,nVerts=0):
    """
    nsurf=1 for first surface,
    nsurf=-1 for second surface
    
    xadd has to be set for second surface (only)
    """

    surfadd=0
    if nsurf == -1:
        surfadd = N2*N1
    
    verts = []
    faces = []
    
    sig = 1
    if rad < 0:
        sig = -1
    rad = np.abs(rad)
    ang = np.arcsin(lrad/rad)

    verts.append(Vector((-xadd,0,0)))
    a = ang/N1
    r = rad*np.sin(a)
    x = rad-np.sqrt(rad**2-r**2)
    for j in range(N2)[::nsurf]:
        b = 2*np.pi*j/N2
        verts.append(Vector((-1.*x*sig*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
        fi1 = nVerts+surfadd
        fi2 = fi1+nsurf*((j+1)%N2+1)
        fi3 = fi1+nsurf*(j+1)
        faces.append([fi1,fi2,fi3])
    for i in range(1,N1):
        a = ang/N1*(i+1)
        r = rad*np.sin(a)
        x = rad-np.sqrt(rad**2-r**2)
        for j in range(N2)[::nsurf]:
            b = 2*np.pi*j/N2
            verts.append(Vector((-1.*x*sig*nsurf-xadd,r*np.sin(b),r*np.cos(b))))
            fi1 = nVerts+surfadd+nsurf*(j+1+i*N2)
            fi2 = nVerts+surfadd+nsurf*((j+1)%N2+1+i*N2)
            fi3 = fi2-nsurf*N2
            fi4 = fi1-nsurf*N2
            faces.append([fi4,fi3,fi2,fi1])
            
    return verts, faces


class OBJECT_OT_add_lens(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}
    
    rad1 = FloatProperty(
           name="Surface 1 Radius",
           default = 12.,
           description="Radius of Curvature of Surface 1",
           )
    rad2 = FloatProperty(
           name="Surface 2 Radius",
           default = 24.,
           description="Radius of Curvature of Surface 2",
           )
    num1 = IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           )
    num2 = IntProperty(
           name="N2",
           default = 64,
           description="Number of angular vertices",
           )
    lensradius = FloatProperty(
           name="Lens Radius",
           default = 3.,
           description="Lens outer radius",
           )
    centerthickness = FloatProperty(
           name="Center Thickness",
           default = 1.,
           description="Center thickness of lens",
           )

    def execute(self, context):

        add_lens(self, context)

        return {'FINISHED'}

class OBJECT_OT_add_mirror(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_mirror"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}
    
    rad = FloatProperty(
           name="Surface 1 Radius",
           default = 12.,
           description="Radius of Curvature of Mirror Surface",
           )
    num1 = IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           )
    num2 = IntProperty(
           name="N2",
           default = 64,
           description="Nubmer of angular vertices",
           )
    mirrorradius = FloatProperty(
           name="Lens Radius",
           default = 3.,
           description="Mirror outer radius",
           )
    centerthickness = FloatProperty(
           name="Thickness",
           default = 1.,
           description="Thickness at thinnest point",
           )

    def execute(self, context):

        add_mirror(self, context)

        return {'FINISHED'}

# Registration

def add_lens_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_lens.bl_idname,
        text="Add Lens",
        icon='PLUGIN')
def add_mirror_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_mirror.bl_idname,
        text="Add Mirror",
        icon='PLUGIN')


# This allows you to right click on a button and link to the manual
def add_lens_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_lens", "editors/3dview/object"),
        )
    return url_manual_prefix, url_manual_mapping
def add_mirror_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_mirror", "editors/3dview/object"),
        )
    return url_manual_prefix, url_manual_mapping


def register():
    bpy.utils.register_class(OBJECT_OT_add_mirror)
    bpy.utils.register_manual_map(add_mirror_manual_map)
    bpy.types.INFO_MT_mesh_add.append(add_mirror_button)
    bpy.utils.register_class(OBJECT_OT_add_lens)
    bpy.utils.register_manual_map(add_lens_manual_map)
    bpy.types.INFO_MT_mesh_add.append(add_lens_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_mirror)
    bpy.utils.unregister_manual_map(add_mirror_manual_map)
    bpy.types.INFO_MT_mesh_add.remove(add_mirror_button)
    bpy.utils.unregister_class(OBJECT_OT_add_lens)
    bpy.utils.unregister_manual_map(add_lens_manual_map)
    bpy.types.INFO_MT_mesh_add.remove(add_lens_button)


if __name__ == "__main__":
    register()
