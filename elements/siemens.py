"""
Copyright 2019-2024, Johannes Hinrichs

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

import bpy
import numpy as np

from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from mathutils import Vector

from .. import surface as sfc
from .. import object_data_add
from .. import utils

class OBJECT_OT_add_siemens(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_siemensstar"
    bl_label = "Siemens Star"
    bl_options = {'REGISTER', 'UNDO'}

    openended : BoolProperty(
            name="Open Ended",
            default=True,
           )
    rad1 : FloatProperty(
           name="Star Radius",
           default = 12.,
           description="Radius of the Siemens star",
           unit = "LENGTH",
           )
    rad2 : FloatProperty(
           name="Outside Radius",
           default = 16.,
           description="Radius of the full disk",
           unit = "LENGTH",
           )
    thickness : FloatProperty(
           name="Thickness",
           default = 2.,
           description="Thickness of the Siemens star",
           unit = "LENGTH",
           )
    num1 : IntProperty(
           name="N",
           default = 32,
           description="Number of line pairs",
           min=4,
           )
    num2 : IntProperty(
           name="Rim Resolution",
           default = 4,
           description="Resolution of the outer rim",
           min=1,
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
        layout.prop(self, 'openended')
        layout.prop(self, 'rad1')
        if not self.openended:
            layout.prop(self, 'rad2')
        layout.prop(self, 'thickness')
        layout.prop(self, 'num1')
        layout.prop(self, 'num2')

    def execute(self, context):
        add_siemensstar(self, context)
        return {'FINISHED'}

def add_siemensstar(self, context):
    verts = []
    edges = []
    faces = []

    #basic parameters
    OE = self.openended
    Nlines = self.num1
    Nseg = self.num2
    Nouter = 2*Nlines*Nseg
    t = self.thickness
    r1 = self.rad1
    r2 = self.rad2

    if not OE:
        if r1 > r2: r2, r1 = r1, r2

    #center vertex
    verts.append(Vector((0,0,0)))

    #Siemens star vertices
    for i in range(2*Nlines):
        alpha = 2*np.pi/2/Nlines*i
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        verts.append(Vector((r1*sa,r1*ca,0)))
        if i%2 == 0:
            for j in range (1,Nseg):
                alpha = 2*np.pi/2/Nlines*(i + j/Nseg)
                ca = np.cos(alpha)
                sa = np.sin(alpha)
                verts.append(Vector((r1*sa,r1*ca,0)))

    nVerts0 = len(verts)
    #outer ring vertices    #if not openended
    if not OE:
        for i in range(Nouter):
            alpha = 2*np.pi/Nouter*i
            ca = np.cos(alpha)
            sa = np.sin(alpha)
            verts.append(Vector((r2*sa,r2*ca,0)))

    #duplicate all vertices for top/bottom
    nVerts = len(verts)
    verts2 = []
    for v in verts:
        verts2.append(Vector((v[0],v[1],t)))
    verts = verts + verts2

    #faces inside star
    for i in range(2*Nlines):
        fi1 = 0
        fi2 = nVerts
        fi3 = 1 + (i+1)//2*Nseg + (i)//2
        fi4 = fi3 + nVerts
        if OE:
            dir = 1 - 2*(i%2)
        else:
            dir = 2*(i%2) - 1
        faces.append([fi1,fi2,fi4,fi3][::dir])
    
    #faces outside star
    for i in range(Nlines):
        for j in range(Nseg):
            f = []
            f.append(i*(Nseg+1) + j + 2)
            f.append(i*(Nseg+1) + j + 1)
            f.append(i*(Nseg+1) + j + 1 + nVerts)
            f.append(i*(Nseg+1) + j + 2 + nVerts)
            if OE:
                faces.append(f)
            else:
                faces.append(f[::-1])
    
    #top/bottom faces
    if OE:
        for i in range(Nlines):
            #bottom
            f = [0]
            for j in range(Nseg+1):
                f.append(i*(Nseg+1) + j + 1)
            faces.append(f)
            #top
            g = [fi + nVerts for fi in f[::-1]]
            faces.append(g)
    else:
        #star pattern
        for i in range(2*Nlines)[1::2]:
            #bottom
            fi1 = 0
            fi2 = 1 + (i+1)//2*Nseg + (i)//2
            fi3 = 1 + ((i+2)//2*Nseg + (i+1)//2)%(nVerts0-1)
            f = [fi1, fi2, fi3]
            faces.append(f)
            #top
            g = [fi + nVerts for fi in f[::-1]]
            faces.append(g)
        #outside ring
        for i in range(2*Nlines):
            #bottom
            f = [nVerts0 + (i*Nseg + j)%Nouter for j in range(Nseg+1)]
            if i%2 == 0:
                fi2 = 1 + (i+1)//2*Nseg + (i)//2
                fi3 = 1 + ((i+2)//2*Nseg + (i+1)//2)%(nVerts0-1)
                f = f + [x for x in np.arange(fi2,fi3+1)[::-1]]
            else:
                fi2 = 1 + (i+1)//2*Nseg + (i)//2
                fi3 = 1 + ((i+2)//2*Nseg + (i+1)//2)%(nVerts0-1)
                f = f + [fi3, fi2]
            faces.append(f)
            #top
            g = [fi + nVerts for fi in f[::-1]]
            faces.append(g)

    #faces of outer ring
    if not OE:
        for i in range(Nouter):
            fi1 = nVerts0 + (i+1)%(Nouter)
            fi2 = nVerts0 + i%(Nouter)
            fi3 = fi2 + nVerts
            fi4 = fi1 + nVerts
            faces.append([fi1,fi2,fi3,fi4])


    mesh = bpy.data.meshes.new(name="New Siemens Star")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

