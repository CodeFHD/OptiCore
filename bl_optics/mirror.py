"""
Copyright 2019-2025, Johannes Hinrichs

This file is part of OptiCore.

OptiCore is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OptiCore is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OptiCore. If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

import bpy
from bpy.types import Operator, ParticleSettingsTextureSlot
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils

class OBJECT_OT_add_mirror(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_mirror"
    bl_label = "Mirror"
    bl_options = {'REGISTER', 'UNDO'}
    
    mtype : EnumProperty(
           name="Surface Shape",
           items = {("parabolic","Parabolic",""),
                    ("spherical","Spherical",""),
                    ("aspheric", "Aspheric", "")},
           default = "parabolic",
           description="Shape of Mirror Surface",
           )
    opos : EnumProperty(
           name="Origin position",
           items = {("FP","Focal Point",""),
                    ("MC","Mirror Center","")},
           default = "FP",
           description="Position of the Mesh Origin w.r.t. optical properties",
           )
    rad : FloatProperty(
           name="Surface Radius",
           default = 12.,
           description="Radius of Curvature of Mirror Surface",
           unit = "LENGTH",
           )
    num1 : IntProperty(
           name="N1",
           default = 32,
           description="Number of radial vertices",
           min=3,
           )
    num2 : IntProperty(
           name="N2",
           default = 64,
           description="Nubmer of angular vertices",
           min=3,
           )
    mirrorradius : FloatProperty(
           name="Mirror Radius",
           default = 3.,
           description="Mirror outer radius",
           unit = "LENGTH",
           )
    centerthickness : FloatProperty(
           name="Back Thickness",
           default = 1.,
           description="Thickness at thinnest point",
           unit = "LENGTH",
           )
    ASPDEG = 3
    k : FloatProperty(
           name="k",
           default = 0.,
           description="Aspheric conical constant",
           )
    A : FloatVectorProperty(
           name="A",
           #default = (0.,0.,0.),
           default = list([0. for i in range(ASPDEG)]),
           description="Aspheric correction coefficients",
           size = ASPDEG,
           )
    theta : FloatProperty(
           name="Offset Angle",
           default = 0.,
           description="Offset angle for off-axis mirror",
           unit = "ROTATION",
           )
    material_name : StringProperty(
            name="Material",
            default="",
           )
    shade_smooth : BoolProperty(
            name="Smooth Shading",
            default=True,
           )
    smooth_type : BoolProperty(
            name="Use Custom Normals",
            default=True,
           )
    cent_hole : BoolProperty(
            name="Central Hole",
            default=False,
           )
    hole_rad : FloatProperty(
           name="HoleRadius",
           default = 0.1,
           description="Radius of Curvature of Mirror Surface",
           min = 0.01,
           unit = "LENGTH",
           )
    display_edit : BoolProperty(
            name="Display Edit Mode",
            default=False,
           )

    def draw(self, context):
        layout = self.layout
        #scene = context.scene
        layout.prop(self, 'mtype')
        if self.mtype == 'parabolic':
            layout.prop(self, 'opos')
        layout.prop(self, 'rad')
        layout.prop(self, 'mirrorradius')
        layout.prop(self, 'centerthickness')
        if self.mtype == 'aspheric':
            layout.prop(self, 'k')
            layout.prop(self, 'A')
        if self.mtype == 'parabolic':
            layout.prop(self, 'theta')
        layout.prop(self, 'num1')
        layout.prop(self, 'num2')
        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        #layout.prop(self, 'material_name')
        layout.prop(self, 'shade_smooth')
        # if self.shade_smooth:
        #     layout.prop(self, 'smooth_type')
        layout.prop(self, 'cent_hole')
        if self.cent_hole:
            layout.prop(self, 'hole_rad')
        layout.prop(self, 'display_edit')

    def execute(self, context):
        add_mirror(self, context)
        return {'FINISHED'}
    
def get_default_paramdict_mirror():
    """
    This functions returns a parameter dictionary filled with some default values.
    This can be used as a template for external calls to add_lens()
    when it is desired to only use a subset of the possible varaibles
    """
    paramdict = {}
    paramdict['rad'] = 12.
    paramdict['num1'] = 32
    paramdict['num2'] = 64
    paramdict['mirrorradius'] = 3.
    paramdict['centerthickness'] = 1.
    paramdict['mtype'] = 'parabolic'
    paramdict['k'] = 0.
    paramdict['A'] = list([0. for i in range(3)])
    paramdict['theta'] = 0.
    paramdict['opos'] = "FP"
    paramdict['cent_hole'] = False
    paramdict['hole_rad'] = 0.1
    paramdict['material_name'] = ""
    paramdict['shade_smooth'] = True
    paramdict['smooth_type'] = True
    paramdict['display_edit'] = False
    paramdict['flipdirection'] = 1
    return paramdict

def add_mirror(self, context, paramdict=None):
    if paramdict is None:
        srad = self.rad
        N1 = self.num1
        N2 = self.num2
        mrad = self.mirrorradius
        CT = self.centerthickness
        surftype = self.mtype
        k = self.k
        A = self.A
        theta = self.theta*180/np.pi
        opos = self.opos
        hrad = 0
        cent_hole = self.cent_hole
        if self.cent_hole:
            hrad = self.hole_rad
            if hrad > 0.99*mrad:
                hrad = 0.99*mrad
        material_name = self.material_name
        shade_smooth = self.shade_smooth
        smooth_type = self.smooth_type
        display_edit = self.display_edit
        flipdirection = False
    else:
        srad = paramdict['rad']
        N1 = paramdict['num1']
        N2 = paramdict['num2']
        mrad = paramdict['mirrorradius']
        CT = paramdict['centerthickness']
        surftype = paramdict['mtype']
        k = paramdict['k']
        A = paramdict['A']
        theta = paramdict['theta']*180/np.pi
        opos = paramdict['opos']
        hrad = 0
        cent_hole = paramdict['cent_hole']
        if cent_hole:
            hrad = paramdict['hole_rad']
            if hrad > 0.99*mrad:
                hrad = 0.99*mrad
        material_name = paramdict['material_name']
        shade_smooth = paramdict['shade_smooth']
        smooth_type = paramdict['smooth_type']
        display_edit = paramdict['display_edit']
        flipdirection = paramdict['flipdirection']
    
    #check surface radius for consistency
    if surftype == 'spherical':
        if abs(srad) < mrad: srad = 0
        # if not utils.check_surface(np.abs(srad), mrad): srad=0

    # no hole for aspheric surface
    if surftype == 'aspheric':
        cent_hole = False

    #compute mirror surface
    if srad == 0: #flat surface case
        verts, faces, normals = sfc.add_flat_surface(mrad, N1, N2, hole=cent_hole, hrad=hrad)
        xadd = 0
        zOA = 0
    else:
        if surftype == 'spherical':
            verts, faces, normals = sfc.add_spherical_surface(-srad, mrad, N1, N2,hole=cent_hole,hrad=hrad)
            xadd = 0
            zOA = 0
            if srad < 0:
                zOA = -np.abs(srad)+np.sqrt(srad**2-mrad**2)
        elif surftype == 'aspheric':
            verts, faces, normals = sfc.add_sagsurface_circular(-srad, k, A, mrad, N1, N2, nVerts=0, surftype='aspheric')
            #verts, faces, normals = sfc.add_aspheric_surface(-srad, k, A, mrad, N1, N2, 1)
            xadd = 0
            zOA = 0
        elif surftype == 'parabolic':
            verts, faces, xadd, zOA, normals = sfc.add_parabolic_surface(-srad, mrad, N1, N2, theta, orig=opos,hole=cent_hole,hrad=hrad)

    # reorder the coordinate system to Blender covnention
    verts = [[t[2], t[0], t[1]] for t in verts]
    normals = [[t[2], t[0], t[1]] for t in normals]

    nVerts = len(verts)
    
    #add rear surface
    verts2, faces2, normals2 = sfc.add_flat_surface(mrad, N1, N2, CT-zOA, xadd, nVerts=nVerts,hole=cent_hole,hrad=hrad)

    # reorder the coordinate system to Blender covnention
    verts2 = [[t[2], t[0], t[1]] for t in verts2]
    normals2 = [[t[2], t[0], t[1]] for t in normals2]

    # flip normals for last surface
    faces2 = [df[::-1] for df in faces2]
    normals2 = [(-t[0], -t[1], -t[2]) for t in normals2]

    nVerts2 = len(verts2)
    
    verts = verts + verts2
    faces = faces + faces2
    normals = normals + normals2

    nVerts_tot = len(verts)
        
    del verts2
    del faces2

    #add side
    for j in range(N2):
        fi1 = nVerts + (j+1)%(N2) + N2*cent_hole
        fi2 = nVerts + j + N2*cent_hole
        fi3 = fi2 - N2*(1 + cent_hole)
        fi4 = fi1 - N2*(1 + cent_hole)
        faces.append([fi4,fi3,fi2,fi1])

    normalsside = sfc.get_ringnormals(N2)
    normalsside = [[t[2], t[0], t[1]] for t in normalsside]

    #fill hole
    normalshole = [] #in case there is no hole
    if cent_hole:
        for j in range(N2):
            fi2 = j
            fi1 = (j+1)%N2
            fi4 = nVerts + (j+1)%N2
            fi3 = nVerts + j
            faces.append([fi4,fi3,fi2,fi1])

        normalshole = sfc.get_ringnormals(N2)
        normalshole = [(-n[2], -n[0], -n[1]) for n in normalshole]

    #create mesh from verts and faces
    mesh = bpy.data.meshes.new(name="Mirror")
    edges = [] # edges are not explicitly created here so we pass an empty list
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)
    
    #assign matierals
    if material_name in bpy.data.materials:
        mat = bpy.data.materials[material_name]
        obj.data.materials.append(mat)

    #apply smooth shading
    if shade_smooth:
        bpy.ops.object.shade_smooth()
        if smooth_type:    #assign custom normals
            bpy.ops.mesh.customdata_custom_splitnormals_clear()
            bpy.ops.mesh.customdata_custom_splitnormals_add()

            cn1, cn2, cn3, cn4 = [], [], [], []
            #mirror surface
            if srad == 0:
                if cent_hole:
                    nloopsface = 4*N2
                else:
                    nloopsface = N2
            else:
                if cent_hole:
                    nloopsface = N2*4*(N1-1)
                else:
                    nloopsface = N2*(3 + 4*(N1-1))
            for i in range(nloopsface):
                vi = mesh.loops[i].vertex_index
                cn1.append(normals[vi])
            #rear surface
            if cent_hole:
                nloopsrear = 4*N2
            else:
                nloopsrear = N2
            for i in range(nloopsrear):
                vi = mesh.loops[i+nloopsface].vertex_index
                cn2.append(normals[vi])
            #side
            nloopsside = 4*N2
            for i in range(nloopsside):
                vi = mesh.loops[i+nloopsface+nloopsrear].vertex_index
                if vi < nVerts:
                    vi = vi - nVerts + N2
                else:
                    vi = (vi - nVerts)%N2
                cn3.append(normalsside[vi])
            #hole
            if cent_hole:
                nloopshole = 4*N2
                for i in range(nloopshole):
                    vi = mesh.loops[i+nloopsface+nloopsrear+nloopsside].vertex_index
                    if vi < nVerts:
                        pass
                    else:
                        vi = vi - nVerts_tot + N2
                    cn4.append(normalshole[vi])

            mesh.normals_split_custom_set(cn1 + cn2 + cn3 + cn4)
            
    if flipdirection == -1:
        obj.rotation_euler = [0,0,np.pi]

    #for testing
    if display_edit:
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')