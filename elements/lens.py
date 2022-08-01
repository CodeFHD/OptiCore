import bpy
import numpy as np

from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add

from .. import surface as sfc
from .. import object_data_add
from .. import utils

class OBJECT_OT_add_lens(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_lens"
    bl_label = "Lens"
    bl_options = {'REGISTER', 'UNDO'}
    
    makedoublet : BoolProperty(
            name="Doublet Lens",
            default=False,
           )
    ltype1 : EnumProperty(
           name="Surface 1 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 1",
           #options={'HIDDEN'},
           )
    ltype2 : EnumProperty(
           name="Surface 2 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 2",
           #options={'HIDDEN'},
           )
    ltype3 : EnumProperty(
           name="Surface 3 Type",
           items = {("spherical","Spherical",""),
                    ("aspheric","Aspheric","")},
           default = "spherical",
           description="Shape of Surface 3",
           #options={'HIDDEN'},
           )
    rad1 : FloatProperty(
           name="Surface 1 Radius",
           default = 12.,
           description="Radius of Curvature of Surface 1",
           unit = "LENGTH",
           )
    rad2 : FloatProperty(
           name="Surface 2 Radius",
           default = -24.,
           description="Radius of Curvature of Surface 2",
           unit = "LENGTH",
           )
    rad3 : FloatProperty(
           name="Surface 3 Radius",
           default = 24.,
           description="Radius of Curvature of Surface 3",
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
           description="Number of angular vertices",
           min=3,
           )
    lensradius : FloatProperty(
           name="Lens Radius",
           default = 3.,
           description="Lens outer radius",
           unit = "LENGTH",
           )
    flangerad1 : FloatProperty(
           name="Flange Width Surface 1",
           default = 0.,
           description="Flange Width Surface 1",
           unit = "LENGTH",
           min = 0,
           )
    flangerad2 : FloatProperty(
           name="Flange Width Surface 2",
           default = 0.,
           description="Flange Width Surface 2",
           unit = "LENGTH",
           min = 0,
           )
    flangerad3 : FloatProperty(
           name="Flange Width Surface 3",
           default = 0.,
           description="Flange Width Surface 3",
           unit = "LENGTH",
           min = 0,
           )
    centerthickness : FloatProperty(
           name="Center Thickness",
           default = 1.,
           description="Center thickness of lens",
           unit = "LENGTH",
           )
    centerthickness2 : FloatProperty(
           name="Center Thickness 2",
           default = 1.,
           description="Center thickness of lens segment 2",
           unit = "LENGTH",
           )
    k : FloatProperty(
           name="k",
           default = 0.,
           description="Aspheric conical constant",
           )
    A : FloatVectorProperty(
           name="A",
           default = (0.,0.,0.),
           description="Aspheric correction coefficients",
           )
    k2 : FloatProperty(
           name="k2",
           default = 0.,
           description="Aspheric conical constant",
           )
    A2 : FloatVectorProperty(
           name="A2",
           default = (0.,0.,0.),
           description="Aspheric correction coefficients",
           )
    k3 : FloatProperty(
           name="k3",
           default = 0.,
           description="Aspheric conical constant",
           )
    A3 : FloatVectorProperty(
           name="A3",
           default = (0.,0.,0.),
           description="Aspheric correction coefficients",
           )
    material_name : StringProperty(
            name="Material",
            default="",
           )
    material_name2 : StringProperty(
            name="Material 2",
            default="",
           )
    material_name3 : StringProperty(
            name="Material 3",
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
    dshape : BoolProperty(
            name="D-Shaped Lens",
            default=False,
           )
    optiverts : BoolProperty(
            name="Alternative radial vertex distribution",
            default=False,
           )
    debugmode : BoolProperty(
            name="Debug Mode",
            default=False,
           )
    ior : FloatProperty(
           name="IOR (Info Only)",
           default = 1.5,
           description="Index of Refraction. Not transferred to material, only for UI information.",
           )
    ior2 : FloatProperty(
           name="IOR (Info Only)",
           default = 1.6,
           description="Index of Refraction of second lens half. Not transferred to material, only for UI information.",
           )
    def get_flen(self):
        return self.flen_intern
    flen_intern : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Focal length of the lens.",
           )
    flen : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Focal length of the lens.",
           get=get_flen,
           )
    addrayfan : BoolProperty(
            name="Add Ray Fan",
            default=False,
           )
    zdet : FloatProperty(
           name="Detector Distance",
           default = 48.0,
           description="Distance to plane where rays are traced to.",
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
        layout.prop(self, 'ltype1')
        layout.prop(self, 'rad1')
        if not self.ltype1=='aspheric': #TMP
            layout.prop(self, 'flangerad1')
        layout.prop(self, 'ltype2')
        layout.prop(self, 'rad2')
        if not self.ltype2=='aspheric': #TMP
            layout.prop(self, 'flangerad2')
        if self.makedoublet:
            layout.prop(self, 'ltype3')
            layout.prop(self, 'rad3')
            if not self.ltype3=='aspheric': #TMP
                layout.prop(self, 'flangerad3')
        layout.prop(self, 'lensradius')
        layout.prop(self, 'centerthickness')
        if self.makedoublet:
            layout.prop(self, 'centerthickness2')
        layout.prop(self, 'num1')
        layout.prop(self, 'num2')
        if self.ltype1=='aspheric':
            layout.prop(self, 'k')
            layout.prop(self, 'A')
        if self.ltype2=='aspheric':
            layout.prop(self, 'k2')
            layout.prop(self, 'A2')
        if self.makedoublet:
            if self.ltype3=='aspheric':
                layout.prop(self, 'k3')
                layout.prop(self, 'A3')
        layout.prop_search(self, "material_name", bpy.data, "materials", icon="NONE")
        if self.makedoublet:
            layout.prop_search(self, "material_name2", bpy.data, "materials", icon="NONE")
            layout.prop_search(self, "material_name3", bpy.data, "materials", icon="NONE")
        layout.prop(self, 'makedoublet')
        layout.prop(self, 'shade_smooth')
        if self.shade_smooth:
            layout.prop(self, 'smooth_type')
        layout.prop(self, 'dshape')
        #layout.prop(self, 'optiverts')
        #layout.prop(self, 'debugmode')
        col2 = layout.column(align=True)
        col2.label(text="Optical Parameters")
        col2.prop(self, 'ior')
        col2.prop(self, 'flen')
        if self.makedoublet:
            col2.prop(self, 'ior2')
        col2.prop(self, 'addrayfan')
        if self.addrayfan:
            col2.prop(self, 'zdet')


    def execute(self, context):
        add_lens(self, context)
        if self.addrayfan:
            utils.trace_rays(self, context)
        return {'FINISHED'}

def add_lens(self, context):
    edges = []
    
    srad1 = self.rad1
    srad2 = -self.rad2
    N1 = self.num1
    N2 = self.num2
    lrad = self.lensradius
    flrad1 = self.flangerad1
    flrad2 = self.flangerad2
    hasfl1 = flrad1 > 0.01*lrad
    hasfl2 = flrad2 > 0.01*lrad
    CT = self.centerthickness
    if self.ltype1 == 'aspheric':
        hasfl1 = 0 #TMP
        k = self.k
        A = self.A
    if self.ltype2 == 'aspheric':
        hasfl2 = 0 #TMP
        k2 = self.k2
        A2 = self.A2
    ssig1 = 1
    if srad1 < 0:
        ssig1 = -1
    ssig2 = 1
    if srad2 < 0:
        ssig2 = -1

    md = self.makedoublet
    if md:
        srad3 = -self.rad3
        flrad3 = self.flangerad3
        hasfl3 = flrad3 > 0.01*lrad
        CT2 = self.centerthickness2
        if self.ltype3 == 'aspheric':
            hasfl3 = 0 #TMP
            k3 = self.k3
            A3 = self.A3
        ssig3 = 1
        if srad3 < 0:
            ssig3 = -1

    #TODO: This may need to move down
    if flrad1 > 0.99*lrad: flrad1 = 0.99*lrad
    if flrad2 > 0.99*lrad: flrad2 = 0.99*lrad
    lrad1 = lrad - flrad1 if hasfl1 else lrad
    lrad2 = lrad - flrad2 if hasfl2 else lrad
    if md:
        if flrad3 > 0.99*lrad: flrad3 = 0.99*lrad
        lrad3 = lrad - flrad3 if hasfl3 else lrad
    
    #check surface radii for consistency
    ##check radius overflow
    if not utils.check_surface(np.abs(srad1), lrad1):
        srad1 = 0
        lrad1 = lrad
        hasfl1 = 0
    if not utils.check_surface(np.abs(srad2), lrad2):
        srad2 = 0
        lrad2 = lrad
        hasfl2 = 0
    if md:
        if not utils.check_surface(np.abs(srad3), lrad3):
            srad3 = 0
            lrad3 = lrad
            hasfl3 = 0
    ##check center thickness
    lsurf1, lsurf2 = 0, 0
    if not srad1 == 0:
        lsurf1 = srad1-ssig1*np.sqrt(srad1**2-lrad1**2)
    if not srad2 == 0:
        lsurf2 = srad2-ssig2*np.sqrt(srad2**2-lrad2**2)
    if (lsurf1 + lsurf2) > CT:
        CT = lsurf1 + lsurf2
    if md:
        lsurf3 = 0
        if not srad3 == 0:
            lsurf3 = srad3-ssig3*np.sqrt(srad3**2-lrad3**2)
        if (-lsurf2 + lsurf3) > CT2: CT2 = -lsurf2 + lsurf3

    #add surface1
    if srad1 == 0: #flat surface case
        verts, faces, normals = sfc.add_flat_surface(lrad1,N1,N2, dshape=self.dshape)
    elif self.ltype1 == 'spherical':
        verts, faces, normals, N1, N2 = sfc.add_spherical_surface(srad1, lrad1, N1, N2, dshape=self.dshape, optiverts=self.optiverts, lrad_ext=lrad)
    elif self.ltype1 == 'aspheric':
        verts, faces, normals = sfc.add_aspheric_surface(srad1, k, A, lrad1, N1, N2, dshape=self.dshape)

    #add flat annulus verts if needed
    #dvert = 
    #dfac = 
    
    nVerts = len(verts)
    nFacs = len(faces)
    
    #add surface2
    if srad2 == 0:#flat surface case
        dvert, dfac, normals2 = sfc.add_flat_surface(lrad2,N1,N2,-1,CT,nVerts=nVerts, dshape=self.dshape)
        #dvert = dvert[::-1]
    elif self.ltype2 == 'spherical':
        dvert, dfac, normals2, N12, N22 = sfc.add_spherical_surface(srad2, lrad2, N1, N2, -1, CT, nVerts=nVerts, dshape=self.dshape, optiverts=self.optiverts, lrad_ext=lrad)
    elif self.ltype2 == 'aspheric':
        dvert, dfac, normals2 = sfc.add_aspheric_surface(srad2, k2, A2, lrad2, N1, N2, -1, CT, nVerts=nVerts, dshape=self.dshape)
    dvert, dfac, normals2 = dvert[::-1], dfac[::-1], normals2[::-1]#for flat surface, there is only one dfac and normals are all same so it doesn't matter to leave it here

    #reverse face normals if doublet lens #TODO: check if needed for convention, in that case add to custom normals
    #if md:
    #    dfac = [f[::-1] for f in dfac]
        
    verts = verts + dvert
    faces = faces + dfac
    normals = normals + normals2

    nVerts2 = len(verts)
    ndVerts2 = len(dvert)
    nFacs2 = len(faces)
    ndFac2 = len(dfac)

    #add side
    for j in range(N2-self.dshape):
        fi1 = nVerts+(j+1)%(N2)
        fi2 = nVerts+j
        fi3 = fi2-N2
        fi4 = fi1-N2
        faces.append([fi1,fi2,fi3,fi4])
    if self.dshape:
        if srad1==0:
            ptss1 = [N2-1,0]
        else:
            ptss1 = [(x+1)*N2 for x in range(N1)[::-1]] + [0] + [x*N2+1 for x in range(N1)]
        if srad2==0:
            ptss2 = [nVerts2-N2, nVerts2-1]
        else:
            ptss2 = [nVerts + x*N2 for x in range(N1)] + [nVerts2-1] + [nVerts + (x+1)*N2-1 for x in range(N1)[::-1]]
        faces.append(ptss1 + ptss2)

    normalsside = sfc.get_ringnormals(N2, dshape=self.dshape)

    nFacSide = N2

    if md:
        #add surface3
        if srad3 == 0:
            #flat surface case
            dvert, dfac, normals3 = sfc.add_flat_surface(lrad3,N1,N2,-1,CT+CT2,nVerts=nVerts2, dshape=self.dshape)
            #dvert = dvert[::-1]
        elif self.ltype3 == 'spherical':
            dvert, dfac, normals3, N12, N22 = sfc.add_spherical_surface(srad3, lrad3, N1, N2, -1, CT+CT2, nVerts=nVerts2, dshape=self.dshape, optiverts=self.optiverts, lrad_ext=lrad)
        elif self.ltype3 == 'aspheric':
            dvert, dfac, normals3 = sfc.add_aspheric_surface(srad3, k3, A3, lrad3, N1, N2, -1, CT+CT2, nVerts=nVerts2, dshape=self.dshape)
        dvert, dfac, normals3 = dvert[::-1], dfac[::-1], normals3[::-1]
        
        verts = verts+dvert
        faces = faces+dfac
        normals = normals + normals3

        nvs1 = N2*N1+1
        if srad1 == 0:
            nvs1 = N2
        nvs2 = N2*N1+1
        if srad2 == 0:
            nvs2 = N2

        nVerts3 = len(verts)
        ndVerts3 = len(dvert)
        nFacs3 = len(faces)
        ndFac3 = len(dfac)

        #add side
        for j in range(N2-self.dshape):
            fi1 = nVerts2+(j+1)%(N2)
            fi2 = nVerts2+j
            fi3 = fi2-nvs2
            fi4 = fi1-nvs2
            faces.append([fi1,fi2,fi3,fi4])
        
        if self.dshape:
            if srad2==0:
                ptss1 = [nvs1 + N2-1, nvs1]
            else:
                ptss1 = [nVerts + (x+1)*N2-1 for x in range(N1)] + [nVerts2-1] + [nVerts + x*N2 for x in range(N1)[::-1]]   #[nvs1 + (x+1)*N2 for x in range(N1)[::-1]] + [nvs1] + [nvs1 + x*N2+1-N2 for x in range(N1)]
            if srad3==0:
                ptss2 = [nVerts3-N2, nVerts3-1]
            else:
                ptss2 = [nVerts2 + x*N2 for x in range(N1)] + [nVerts3-1] + [nVerts2 + (x+1)*N2-1 for x in range(N1)[::-1]]
            faces.append(ptss1 + ptss2)
    
    #create mesh from verts and faces
    del dvert
    del dfac
    mesh = bpy.data.meshes.new(name="New Lens")
    mesh.from_pydata(verts, edges, faces)
    obj = object_data_add(context, mesh, operator=self)

    #assign material(s)
    if self.material_name in bpy.data.materials:
        mat = bpy.data.materials[self.material_name]
        obj.data.materials.append(mat)
    if md:
        cond2 = self.material_name2 in bpy.data.materials
        cond3 = self.material_name3 in bpy.data.materials
    if md and cond2 and cond3:
        if self.material_name2 in bpy.data.materials:
            mat2 = bpy.data.materials[self.material_name2]
            obj.data.materials.append(mat2)
        if self.material_name3 in bpy.data.materials:
            mat3 = bpy.data.materials[self.material_name3]
            obj.data.materials.append(mat3)
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')
        sel_mode = bpy.context.tool_settings.mesh_select_mode
        bpy.context.tool_settings.mesh_select_mode = [False, False, True]
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        #rear surface (mat3)
        for i in range(nFacSide):
            mesh.polygons[i + nFacs3].select=True
        for i in range(ndFac3):
            mesh.polygons[i + nFacs2+ nFacSide].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        obj.active_material_index = 2
        bpy.ops.object.material_slot_assign()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        #middle surface (mat2)
        for i in range(ndFac2):
            mesh.polygons[i + nFacs].select=True
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.context.tool_settings.mesh_select_mode = sel_mode
        obj.active_material_index = 1
        bpy.ops.object.material_slot_assign()
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

    #apply smooth shading
    if self.shade_smooth:
        bpy.ops.object.shade_smooth()
        mesh.use_auto_smooth = True
        if self.smooth_type: #assign custom normals
            bpy.ops.mesh.customdata_custom_splitnormals_clear()
            bpy.ops.mesh.customdata_custom_splitnormals_add()
            cn1, cn2, cn3, cn4, cn5 = [], [], [], [], []
            #surf1
            if srad1 == 0:
                nloops1 = N2
                nloops1_fl = N2
            else:
                nloops1 = (N2-self.dshape)*(3 + 4*(N1-1))
                nloops1_fl = nloops1 - 2*hasfl1*(3 + 4*(N1-1)+1)
            for i in range(nloops1_fl):
                vi = mesh.loops[i].vertex_index
                cn1.append(normals[vi])
            if hasfl1:
                for i in range(2*(3 + 4*(N1-1))+2):
                    vi = mesh.loops[nloops1_fl + i].vertex_index
                    cn1.append(normals[nVerts-1])
            #surf2
            if srad2 == 0:
                nloops2 = N2
                d_fl = 0
                nloops2_fl = N2
            else:
                nloops2 = (N2-self.dshape)*(3 + 4*(N1-1))
                d_fl = 2*hasfl2*(3 + 4*(N1-1)+1)
                nloops2_fl = nloops2 - 2*hasfl2*(3 + 4*(N1-1)+1)
            if hasfl2:
                for i in range(d_fl):
                    vi = mesh.loops[nloops2_fl + i + nloops1].vertex_index
                    cn2.append(normals[nVerts2-1])
            for i in range(d_fl, nloops2):
                vi = mesh.loops[i + nloops1].vertex_index
                cn2.append(normals[vi])
            #edge12
            nloops3 = 4*(N2-self.dshape)
            for i in range(nloops3):
                vi = mesh.loops[i + nloops1 + nloops2].vertex_index
                if vi < nVerts:
                    vi = vi - nVerts + N2
                else:
                    vi = vi - nVerts
                cn3.append(normalsside[vi])
            #Dface
            if self.dshape:
                if srad1 == 0:
                    nloopsd1 = 1
                else:
                    nloopsd1 = 2*N1
                if srad2 == 0:
                    nloopsd2 = 1
                else:
                    nloopsd2 = 2*N1
                for i in range(nloopsd1 + nloopsd2 + 2):
                    cn3.append((0,-1,0))
                nloops3 += nloopsd1 + nloopsd2 + 2
            if md:
                #surf3
                if srad3 == 0:
                    nloops4 = N2
                    d_fl = 0
                    nloops4_fl = N2
                else:
                    nloops4 = (N2-self.dshape)*(3 + 4*(N1-1))
                    d_fl = 2*hasfl3*(3 + 4*(N1-1)+1)
                    nloops4_fl = nloops4 - 2*hasfl3*(3 + 4*(N1-1)+1)
                if hasfl3:
                    for i in range(d_fl):
                        vi = mesh.loops[nloops4_fl + i + nloops1 + nloops2 + nloops3].vertex_index
                        cn4.append(normals[nVerts3-1])
                for i in range(d_fl, nloops4):
                    vi = mesh.loops[i + nloops1 + nloops2 + nloops3].vertex_index
                    cn4.append(normals[vi])
                """
                #surf3
                if srad3 == 0:
                    nloops4 = N2
                else:
                    nloops4 = (N2-self.dshape)*(3 + 4*(N1-1))
                for i in range(nloops4):
                    vi = mesh.loops[i + nloops1 + nloops2 + nloops3].vertex_index
                    cn4.append(normals[vi])
                """
                #edge23
                nloops5 = 4*(N2-self.dshape)
                for i in range(nloops5):
                    vi = mesh.loops[i + nloops1 + nloops2 + nloops3 + nloops4].vertex_index
                    if vi < nVerts2:
                        vi = vi - nVerts
                    else:
                        vi = vi - nVerts2
                    cn5.append(normalsside[vi])
                #Dface2
                if self.dshape:
                    if srad3 == 0:
                        nloopsd3 = 1
                    else:
                        nloopsd3 = 2*N1
                    for i in range(nloopsd2 + nloopsd3 + 2):
                        cn5.append((0,-1,0))
            mesh.normals_split_custom_set(cn1 + cn2 + cn3 + cn4 + cn5)

    #compute optical parameters
    if not md:
        self.flen_intern = utils.lens_math.f_lensmaker(self.rad1, self.rad2, self.ior, self.centerthickness)
            
    #for testing
    if self.debugmode:
        mesh.calc_normals_split()
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.mesh.select_all(action='DESELECT')