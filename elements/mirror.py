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
    
    #check surface radius for consistency
    if not utils.check_surface(np.abs(srad), mrad): srad=0

    #compute mirror surface
    if srad == 0: #flat surface case
        verts, faces = sfc.add_flat_surface(mrad,N1,N2)
    else:
        verts, faces = sfc.add_spherical_surface(-1.*srad, mrad, N1, N2)
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
    dvert, dfac = sfc.add_flat_surface(mrad,N1,N2,-1,CT,nVerts)
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