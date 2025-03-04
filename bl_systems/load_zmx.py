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

import time
import numpy as np

from mathutils import Vector

from ..utils import debugprint

from ..raytrace.trace_sequential import exec_trace
from ..bl_raytrace.trace_scene import trace_to_scene
from ..raytrace import rayfan
from ..bl_optics import add_lens, add_mirror, get_default_paramdict_lens, get_default_paramdict_mirror, add_circular_aperture, add_sensor
from ..fileparser.zmx import load_from_zmx
from ..bl_materials import glass_from_Element_cycles, clear_all_materials, add_blackoutmaterial_cycles, add_diffusematerial_cycles

import bpy
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add


def reset_loadzmx_faults():
    bpy.context.window_manager.operator_properties_last("mesh.load_zmx").num1 = 8
    bpy.context.window_manager.operator_properties_last("mesh.load_zmx").num2 = 16
    bpy.context.window_manager.operator_properties_last("mesh.load_zmx").nrays = 11
    bpy.context.window_manager.operator_properties_last("mesh.load_zmx").addrayfan = False
    bpy.context.window_manager.operator_properties_last("mesh.load_zmx").tracetoscene = False

""" Menu registered function """

class OBJECT_OT_load_zmx(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.load_zmx"
    bl_label = "Load from .zmx"
    bl_options = {'REGISTER', 'UNDO'}
    
    path: bpy.props.StringProperty(
        name = "File Path",
        description="Choose a file:",
        subtype='FILE_PATH')
    
    wl: bpy.props.FloatProperty(
        name = "Wavelength",
        description="Wavelength in micrometers",
        default = 0.5876,
        min = 0.35,
        max = 2.0,
        precision=3)
    
    dshape : BoolProperty(
            name="Cross-section Model",
            default=False,
           )
    
    num1 : IntProperty(
           name="N1",
           default = 8,
           description="Number of radial vertices",
           min=3,
           )
    num2 : IntProperty(
           name="N2",
           default = 16,
           description="Number of angular vertices",
           min=3,
           )
    def get_flen(self):
        return self.flen_intern
    flen_intern : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Paraxial EFL of the lens",
           )
    flen : FloatProperty(
           name="Focal Length",
           default = 0.,
           description="Paraxial EFL of the lens",
           get=get_flen,
           )
    def get_rmsspotsize(self):
        return self.rmsspotsize_intern
    rmsspotsize_intern : FloatProperty(
           name="RMS spot size",
           default = 0.,
           description="RMS spot size of the ray trace in the sensor plane",
           )
    rmsspotsize : FloatProperty(
           name="RMS spot size",
           default = 0.,
           description="RMS spot size of the ray trace in the sensor plane",
           get=get_rmsspotsize,
           )
    addlenses : BoolProperty(
            name="Add Lenses",
            default=True,
           )
    addaperture : BoolProperty(
            name="Add Aperture",
            default=True,
           )
    addrayfan : BoolProperty(
            name="Add Ray Fan",
            default=False,
           )
    addsensor : BoolProperty(
            name="Add Sensor",
            default=True,
           )
    thicksensor : BoolProperty(
            name="Thick",
            default=True,
           )
    sensorfactor : FloatProperty(
           name="Sensor size factor",
           default = 1.0,
           description="Factor to adjust sensor size",
           min = 0.01,
           max = 3.0,
           )
    zdet : FloatProperty(
           name="Detector Offset",
           default = 0.0,
           description="Displacement of the sensor plane from zmx file",
           )
    nrays : IntProperty(
           name="Nrays",
           default = 11,
           description="Number of rays for ray fan",
           min=3,
           )
    fantype : EnumProperty(
           name="Fan Type",
           items = {("2D","2D",""),
                    ("3D_tri","3D tris",""),
                    ("3D_tri_finite","3D tris - finite",""),
                    ("3D_random","3D random",""),
                    ("3D_random_finite","3D random - finite",""),
                    ("2D_finite","2D - finite",""),
                    ("3D_rings", "3D rings", ""),
                    ("3D_rings_finite", "3D rings - finite", ""),},
           default = "2D",
           description="Ray Fan Type",
           #options={'HIDDEN'},
           )
    fandist : FloatProperty(
           name="Fan Origin",
           default = 20.,
           description="Distance where ray fan originates.",
           )
    fanangle1 : FloatProperty(
           name="Fan FoV Angle",
           default = 0.,
           description="Field-of-View angle of the ray fan.",
           min = -np.pi/2,
           max = np.pi/2,
           unit = "ROTATION",
           )
    fanangle2 : FloatProperty(
           name="Fan Azimuth Angle",
           default = np.pi/2,
           description="Azimuth angle of the ray fan.",
           min = 0,
           max = 2*np.pi,
           unit = "ROTATION",
           )
    fanangle3 : FloatProperty(
           name="Fan Roll Angle",
           default = np.pi/2.,
           description="Roll angle of the ray fan.",
           min = 0.,
           max = np.pi,
           unit = "ROTATION",
           )
    fanangle_additional : StringProperty(
           name="Additional",
           default = "",
           description="Additional fan angles. Enter float values separated by semicolon",
           )
    fandiam : FloatProperty(
           name="Fan diameter factor",
           default = 1.0,
           description="Factor to adjust stop size",
           min = 0.1,
           max = 3.0,
           )
    aperturediam : FloatProperty(
           name="Stop diameter factor",
           default = 1.0,
           description="Factor to adjust stop size",
           min = 0.01,
           max = 3.0,
           )
    onlycompleterays : BoolProperty(
            name="Display only uninterrupted rays",
            default=False,
           )
    tracetoscene : BoolProperty(
            name="Trace Ray Fan to Scene",
            default=False,
           )
    ghost_order : StringProperty(
            name="Ghost sequence",
            description="Specify two indices, seprated by a comma, between which surfaces a ghost reflection shall be traced. First surface has index 1.",
            default="",
           )
    autofocus : BoolProperty(
            name="RMS-Focus",
            description='Attempt to autofocus by ray tracing and evaluating the RMS spot size',
            default=False,
           )
    excludedetector : BoolProperty(
            name="Exclude Detector",
            default=False,
           )
    display_edit : BoolProperty(
        name="Display Edit Mode",
        default=False,
        )
    split_cemented : BoolProperty(
        name="Split models of cemented lenses",
        default=False,
        )
    mat_refract_only : BoolProperty(
        name="Material: Set refractive only",
        default=False,
        )

    def draw(self, context):
        scene = context.scene
        row = self.layout
        col = row.column(align=True)
        col.label(text="Import File")
        col.prop(self, 'path', text="Select File")
        col.label(text="Ray Tracing")
        col.prop(self, 'flen')
        col.prop(self, 'rmsspotsize')
        col.prop(self, 'wl', text="Wavelength (um)")
        col.prop(self, 'zdet')
        col.prop(self, 'nrays')
        col.prop(self, 'fantype')
        col.prop(self, 'fandist')
        col.prop(self, 'fanangle1')
        col.prop(self, 'fanangle2')
        col.prop(self, 'fanangle3')
        col.prop(self, 'fanangle_additional')
        col.prop(self, 'fandiam')
        col.prop(self, 'aperturediam')
        col.prop(self, 'excludedetector')
        col.prop(self, 'addlenses')
        col.prop(self, 'addaperture')
        col.prop(self, 'addrayfan')
        row = col.row()
        row.prop(self, 'addsensor')
        row.prop(self, 'thicksensor')
        col.prop(self, 'sensorfactor')
        col.prop(self, 'tracetoscene')
        col.prop(self, 'ghost_order')
        col.prop(self, 'autofocus')
        col.prop(self, 'onlycompleterays')
        col.label(text="Mechanical Parameters")
        col.prop(self, 'dshape')
        col.prop(self, 'display_edit')
        col.prop(self, 'split_cemented')
        col.prop(self, 'num1')
        col.prop(self, 'num2')
        col.prop(self, 'mat_refract_only')

    
    def execute(self, context):
        t0 = time.perf_counter()
        if self.path == '':
            return {'FINISHED'}
        if not str(self.path.lower()).endswith('.zmx'):
            return {'FINISHED'}
        fname = self.path.split("\\")[-1]
        print()
        print(f'Now importing {fname}')
        # load the lens
        lens = load_from_zmx(bpy.path.abspath(self.path))
        lens.apertures[0]['radius'] = lens.apertures[0]['radius']*self.aperturediam
        created_objects = [] # list of the individual objects that were created
        t1 = time.perf_counter()

        using_cycles = context.scene.render.engine == 'CYCLES'

        # remove orphaned OC_ materials
        if using_cycles:
            clear_all_materials()
        
        # create meshes from the elements
        if self.addlenses:
            CT_sum = 0
            radius_outer = 0
            for ele, dz in lens.elements:
                i0 = 0 # starting index for surface, relevant for split_cemented
                """
                TODO: Possibly to be foreseen in this implementation:
                In case of a multi-bounce system, i.e. where the same surface is used multiple times,
                implement a logic that it is not necessary to create the same element twice.
                """
                num_surfaces = len(ele.data['radius'])
                if num_surfaces > 4 and not self.split_cemented:
                    # TODO: Support more than triplet groups? Realistically needed?
                    continue
                dz_this = CT_sum + dz
                CT_sum = CT_sum + dz + sum(ele.data['CT'][:-1])
                if ele.ismirror:
                    paramdict = get_default_paramdict_mirror()
                    paramdict['rad'] = -1*ele.data['radius'][0]*ele.direction
                    paramdict['num1'] = self.num1
                    paramdict['num2'] = self.num2
                    paramdict['mirrorradius'] = ele.data['rCA']
                    paramdict['centerthickness'] = 1.
                    paramdict['mtype'] = 'aspheric'
                    paramdict['k'] = ele.data['asph'][0]
                    paramdict['A'] = ele.data['asph'][1:]
                    paramdict['theta'] = 0.
                    #paramdict['opos'] = "FP"
                    paramdict['cent_hole'] = False
                    paramdict['hole_rad'] = 0.1
                    paramdict['material_name'] = ""
                    paramdict['shade_smooth'] = True
                    paramdict['smooth_type'] = True
                    paramdict['display_edit'] = False
                    paramdict['flipdirection'] = ele.direction
                    add_mirror(self, context, paramdict=paramdict)
                    bpy.ops.transform.translate(value=(-dz_this, 0, 0))
                    obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                    created_objects.append(obj_name)
                else:
                    if self.split_cemented:
                        n_loop = 2
                    else:
                        n_loop = num_surfaces
                    # get materials
                    if using_cycles:
                        materials_bulk, materials_interface = glass_from_Element_cycles(ele, self.wl, self.mat_refract_only)
                        material_edge = add_blackoutmaterial_cycles()
                        material_dface = add_diffusematerial_cycles(viewportcolor=[1, 0, 0, 1])
                    while True:
                        if self.split_cemented:
                            i1 = i0 + 2
                        else:
                            i1 = num_surfaces+1
                        paramdict = get_default_paramdict_lens()
                        for i in range(n_loop):
                            ry = ele.data['radius'][i0+i]
                            ky = ele.data['asph'][i0+i][0]
                            Ay = ele.data['asph'][i0+i][1:]
                            rz = ele.data['radius2'][i0+i]
                            ltype = ele.data['lenstype'][i0+i]
                            paramdict[f'ltype{i+1}'] = ltype
                            paramdict[f'RY{i+1}'] = ry
                            paramdict[f'k{i+1}'] = ky
                            paramdict[f'A{i+1}'] = Ay
                            paramdict[f'RZ{i+1}'] = rz
                            paramdict[f'surfrot{i+1}'] = ele.data['surf_rotation'][i0+i]
                            if ele.data['rCA_short'][i0+i] < max(ele.data['rCA'][i0:i1]):
                                paramdict[f'flangerad{i+1}'] = max(ele.data['rCA'][i0:i1]) - ele.data['rCA_short'][i0+i]
                            if using_cycles:
                                if i == 0:
                                    paramdict[f'material_name{i+1}'] = materials_bulk[i]
                                elif i == n_loop-1:
                                    paramdict[f'material_name{i+1}'] = materials_bulk[i-1]
                                else:
                                    paramdict[f'material_name{i+1}'] = materials_interface[i-1]
                        for i in range(n_loop-1):
                            paramdict[f'centerthickness{i+1}'] = ele.data['CT'][i0+i]
                        
                        if using_cycles:
                            paramdict[f'material_edge'] = material_edge# 'OC_LensEdge_cycles'
                            paramdict[f'material_dface'] = material_dface# 'OC_LensDface_cycles'
                        paramdict['lensradius'] = max(ele.data['rCA'][i0:i1])
                        paramdict['makedoublet'] = str(n_loop-1)
                        # parameters completely independent of split_cemented
                        paramdict['num1'] = self.num1
                        paramdict['num2'] = self.num2
                        paramdict['dshape'] = self.dshape
                        if ele.outline_shape == 'square':
                            paramdict['squarelens'] = True
                        add_lens(self, context, paramdict=paramdict)
                        bpy.ops.transform.translate(value=(-dz_this, 0, 0))
                        obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                        created_objects.append(obj_name)

                        if i1 >= num_surfaces:
                            #last surface was included
                            # check for >= and not == to be safe against accidental double increments creating infinite loops
                            break
                        else:
                            dz_this = dz_this + ele.data['CT'][i0]
                            i0 = i0 + 1
                radius_outer = max(radius_outer, max(ele.data['rCA']))
                
        t2 = time.perf_counter()

        # apply sensor offset
        lens.detector['distance'] = lens.detector['distance'] + self.zdet
        # build the lens, raytrace, and get the rays for plotting
        lens.build(self.wl)
        self.flen_intern = lens.EFL_paraxial()
        # print('PARAXIAL BFL:', lens.BFL_paraxial())
        
        t3 = time.perf_counter()

        # create meshes for apertures
        # radius_outer = max([r for r in ele.data['rCA'] for ele, dz in lens.elements])
        # create a viewport black shaded material for the aperture
        if self.addaperture:
            for i, ap in lens.apertures.items():
                radius_inner = lens.apertures[0]['radius']
                idx_ap = max(ap['idx_surface'], 1)
                radius_outer = lens.data['lrad'][idx_ap]
                radius_outer = max(radius_inner*1.2, radius_outer)
                add_circular_aperture(self, context, radius_inner, radius_outer, self.num2, dshape = self.dshape)
                bpy.ops.transform.translate(value=(-ap['z_ap'], 0, 0))
                # obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                # created_objects.append(obj_name)
            
        t4 = time.perf_counter()
        t41_sum = 0
        t42_sum = 0
            
        if self.addrayfan:
            t40 = time.perf_counter()
            # try to convert the angles and clamp between -90 and 90 degrees
            fanangles = [self.fanangle1]
            try:
                if not self.fanangle_additional in ["", "None", "none", "NONE"]:
                    fanangles_add = [max(min(float(a), 89.99), -89.99) for a in self.fanangle_additional.split(";")]
                    fanangles = fanangles_add + fanangles
            except:
                print("FANANGLE ERROR")

            for fanangle in fanangles:
                # set up the rays
                initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist, self.fanangle1, self.fanangle2, self.fanangle3] 
                rays = rayfan.RayFan(self.fantype, initparams, store_history=True)
                # try to parse ghost_oder
                ghost_order = self.ghost_order
                nos = lens.num_optical_surfaces
                custom_surflist_valid = False
                if not ghost_order == '':
                    try:
                        surflist = [int(_) for _ in ghost_order.split(',')]
                        if np.any(np.array(surflist) > nos) or np.any(np.array(surflist) < 1): 
                            print("Warning: Custom ghost list invalid. Tracing default!")
                        elif surflist[0] == surflist[1]:
                            pass 
                        else:
                            s0 = min(surflist[:2])
                            s1 = max(surflist[:2])
                            seq1 = [_ for _ in range(1, s1 + 1)] 
                            seq2 = [_ for _ in range(s0 + 1, s1)][::-1]
                            seq3 = [_ for _ in range(s0, nos + 1)] 
                            surflist = seq1 + seq2 + seq3
                            custom_surflist_valid = True
                    except:
                        print("Warning: Custom ghost list invalid. Tracing default!")
                if not custom_surflist_valid:
                    surflist = [i for i in range(1, nos + 1)] # standard surface list for non-ghost trace
                lens.surf_sequence = surflist
                debugprint('surflist', surflist)
                trace_detector = not self.tracetoscene and not self.excludedetector
                rays = exec_trace(lens, rays, surflist, trace_detector=trace_detector)
                if trace_detector and self.autofocus:
                    retval = rays.autofocus(EFL=self.flen_intern)
                    if retval is not None:
                        offset_min, P_new = retval
                        rays.O_history[len(rays.O_history)-1] = P_new
                    else:
                        print('AUTOFOCUS FAILED')
                self.rmsspotsize_intern = rays.calc_rms_spotsize()
        
                # if trace to scene, perform the last trace
                if self.tracetoscene:
                    rays = trace_to_scene(context, rays)
                
                t41 = time.perf_counter()
                t41_sum = t41_sum + t41 - t40
        
                # create the lines
                verts = []
                edges = []
                faces = []
                nVerts = 0
                ocr = self.onlycompleterays
                if ocr:
                    Ohist = np.array([ooo for ooo in rays.O_history.values()]) # need to convert the dict for the next function to work
                    # explanation: indices where the sum over Points is not nan.
                    # Using sum and not nansum makes this nan if any of the points are nan
                    idx_rays = np.where(np.isnan(np.sum(Ohist, axis=0).sum(axis=1)))[0]
                else:
                    idx_rays = range(rays.O_history[0].shape[0])
                #print('IDXRAYS', idx_rays)
                for j in range(len(rays.O_history)-1):
                    i_valid = -1
                    i_solo = 0
                    for i in range(rays.O_history[0].shape[0]):
                        if i in idx_rays and ocr: 
                            #print('discarding ray', j, i)
                            continue
                        # change axis convention from raytrace to Blender orientation
                        o1 = np.array(rays.O_history[j][i])
                        o1[[0,1,2]] = o1[[2,0,1]]
                        o1[0] = -1*o1[0]
                        o2 = np.array(rays.O_history[j+1][i])
                        o2[[0,1,2]] = o2[[2,0,1]]
                        o2[0] = -1*o2[0]
                        anynan = np.any(np.isnan(o1)) or np.any(np.isnan(o2))
                        if not anynan:
                            i_valid = i_valid + 1
                            verts.append(Vector(o1))
                            verts.append(Vector(o2))
                            edges.append([nVerts + 2*i_valid + i_solo, nVerts + 2*i_valid + 1 + i_solo])
                        elif j==0:
                            i_solo = i_solo + 1
                            verts.append(Vector(o1))
                    nVerts = len(verts)
                #print('RAY LOOP DONE')
                mesh = bpy.data.meshes.new(name="Rayfan")
                mesh.from_pydata(verts, edges, faces)
                obj = object_data_add(context, mesh, operator=self)
                t42 = time.perf_counter()
                t42_sum = t42_sum + t42 - t41
            
        #t5 = time.perf_counter()
        """
        debugprint()
        print('Execution Times:')
        print(f'parsing .zmx-file: {(t1-t0):.3f}')
        print(f'Create Element-Meshes: {(t2-t1):.3f}')
        print(f'Build lens: {(t3-t2):.3f}')
        print(f'Create Aperture-Meshes: {(t4-t3):.3f}')
        if self.addrayfan:
            print(f'Raytracing: {(t41_sum):.3f}')
            print(f'Adding ray-fan-mesh: {(t42_sum):.3f}')
        debugprint()
        """

        if self.addsensor:
            lx = lens.detector['sizex']/2*self.sensorfactor
            ly = lens.detector['sizey']/2*self.sensorfactor
            zsensor = lens.data['CT_sum'][-1] + lens.detector['distance']
            add_sensor(self, context, lx, ly, zsensor, thicksensor=self.thicksensor)
        
        # Select all just created objects
        bpy.ops.object.select_all(action='DESELECT')
        for o in bpy.data.objects:
            if o.name in created_objects:
                o.select_set(True)
        
        if self.display_edit:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='SELECT')

        print('Lens import finished. For transfer to Blender materials,\nuse the following refractive indices')
        print('Surface no.\tn_Blender')
        for i in range(1, lens.num_surfaces):
            n1 = lens.data['n'][i-1]
            n2 = lens.data['n'][i]
            if n2 > 1.1:
                nratio = n2/n1
            else:
                nratio = n1/n2
            print(f'{i}\t\t{nratio:.6f}')
        
        return {'FINISHED'}
