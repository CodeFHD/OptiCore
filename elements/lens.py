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
        verts, faces = sfc.add_flat_surface(lrad,N1,N2)
    else:
        verts, faces = sfc.add_spherical_surface(srad1, lrad, N1, N2)
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
        dvert, dfac = sfc.add_flat_surface(lrad,N1,N2,-1,CT,nVerts=nVerts)
        dvert = dvert[::-1]
    else:
        dvert, dfac = sfc.add_spherical_surface(srad2, lrad, N1, N2,-1, CT, nVerts=nVerts)
        dvert, dfac = dvert[::-1], dfac[::-1]
        
    verts = verts+dvert
    faces = faces+dfac
        
    del dvert
    del dfac

    mesh = bpy.data.meshes.new(name="New Lens")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)
    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if self.shade_smooth:
        bpy.ops.object.shade_smooth()
    if self.split_edge:
        bpy.ops.object.modifier_add(type='EDGE_SPLIT')