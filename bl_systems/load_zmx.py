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
from ..bl_optics import add_lens, add_mirror, get_default_paramdict_lens, get_default_paramdict_mirror, add_sensor
from ..bl_optomech.lens_aperture import add_circular_aperture
from ..bl_optomech.lens_housing import add_lenshousing_simple
from ..fileparser.zmx import load_from_zmx
from ..bl_materials import *
from ..bl_scenes.luxcore_lights import add_laser_array

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
    
    display_option : EnumProperty(
        name="Settings",
        items = {("default","Import",""),
                 ("luxcore","LuxCore Features","")},
        default = "default",
        description="Select which settings are displayed in this window.",
           )
    
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
    addhousing : BoolProperty(
            name="Simple Housing",
            default=False,
            description="Adds a very simple stray light housing around the lens, "
            "leaving 99 percent of the lens clear aperture. be advised that it may leak light, "
            "can overlap for certain geometries and manual verification of its performance is advised. "
            "Solidification via modifier might be suitable.",
           )
    housingtype : BoolProperty(
            name="Tight",
            default=False,
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
    sensor_offset : FloatProperty(
           name="Sensor Offset",
           default = 0.0,
           description="Displacement of the sensor plane from zmx file",
           )
    addcamera : BoolProperty(
            name="Add Camera",
            default=True,
            description="Add a camera looking at the sensor. Only applies if a sensor is added",
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
        default=True,
        )
    origin_position : EnumProperty(
        name="Origin",
        items = {("first_vertex","First Lens",""),
                ("sensor","Sensor",""),},
        default = "first_vertex",
        description="Determines which component of the lens is placed at the location of the 3D cursor",
        #options={'HIDDEN'},
        )
    
    # LuxCore Features
    addlasers : BoolProperty(
        name="Add Laser Lights",
        default=False,
        description='Adds a series of laser lights aimed at the front lens to quickly visualize Bokeh.'
        )
    HFOV : FloatProperty(
           name="H-FoV",
           default = 0.0,
           description="Horizontal Field of View spread",
           min = 0,
           max = 90,
           )
    VFOV : FloatProperty(
           name="V-FoV",
           default = 0.0,
           description="Vertical Field of View spread",
           min = 0,
           max = 90,
           )
    NHFOV : IntProperty(
           name="Num-H",
           default = 4,
           description="Number of Horizontal points",
           min = 1,
           soft_max = 10,
           )
    NVFOV : IntProperty(
           name="Num-V",
           default = 4,
           description="Number of Vertical points",
           min = 1,
           soft_max = 10,
           )

    def draw(self, context):
        disp = self.display_option
        scene = context.scene
        row = self.layout
        col = row.column(align=True)
        col.label(text="Import File")
        col.prop(self, 'path', text="Select File")
        if disp == 'default':
            col.label(text="Ray Tracing")
            col.prop(self, 'flen')
            col.prop(self, 'rmsspotsize')
            col.prop(self, 'wl', text="Wavelength (um)")
            col.prop(self, 'nrays')
            col.prop(self, 'fantype')
            col.prop(self, 'fandist')
            col.prop(self, 'fanangle1')
            col.prop(self, 'fanangle2')
            col.prop(self, 'fanangle3')
            col.prop(self, 'fanangle_additional')
            col.prop(self, 'fandiam')
            col.prop(self, 'aperturediam')
            col.prop(self, 'addlenses')
            col.prop(self, 'addaperture')
            row = col.row()
            row.prop(self, 'addhousing')
            row.prop(self, 'housingtype')
            row = col.row()
            row.prop(self, 'addrayfan')
            row.prop(self, 'excludedetector')
            col.prop(self, 'tracetoscene')
            row = col.row()
            row.prop(self, 'addsensor')
            row.prop(self, 'thicksensor')
            col.prop(self, 'sensor_offset')
            col.prop(self, 'sensorfactor')
            col.prop(self, 'addcamera')
            col.prop(self, 'ghost_order')
            col.prop(self, 'autofocus')
            col.prop(self, 'onlycompleterays')
            col.label(text="Mechanical Parameters")
            col.prop(self, 'dshape')
            col.prop(self, 'display_edit')
            col.prop(self, 'split_cemented')
            col.prop(self, 'num1')
            col.prop(self, 'num2')

        elif disp == 'luxcore':
            col.prop(self, 'addlasers')
            row = col.row()
            row.prop(self, 'HFOV')
            col.prop(self, 'NHFOV')
            row = col.row()
            row.prop(self, 'VFOV')
            col.prop(self, 'NVFOV')

        # global config, always displayed
        col.prop(self, 'mat_refract_only')
        col.prop(self, 'origin_position')
        col.prop(self, 'display_option')

    
    def execute(self, context):
        t0 = time.perf_counter()
        if self.path == '':
            return {'FINISHED'}
        if not str(self.path.lower()).endswith('.zmx'):
            return {'FINISHED'}
        fname = self.path.split("\\")[-1]
        print(f'[OC] Now importing .zmx-file: {fname}')
        # load the lens
        lens = load_from_zmx(bpy.path.abspath(self.path))
        lens.apertures[0]['radius'] = lens.apertures[0]['radius']*self.aperturediam
        lens.detector['npixx'] = int(lens.detector['npixx']*self.sensorfactor)
        lens.detector['npixy'] = int(lens.detector['npixy']*self.sensorfactor)
        lens.detector['sizex'] = lens.detector['npixx']*lens.detector['pixelpitch']
        lens.detector['sizey'] = lens.detector['npixy']*lens.detector['pixelpitch']

        created_lenses = [] # list of (object names) of the individual lenses
        created_rayfans = [] # list of (object names) of the ray fans
        created_other = [] # list of (object names) for other components: aperture, housing, sensor, camera
        created_lights = []
        t1 = time.perf_counter()

        using_cycles = context.scene.render.engine == 'CYCLES'
        using_luxcore = context.scene.render.engine == 'LUXCORE'

        # remove orphaned OC_ materials
        if using_cycles or using_luxcore:
            clear_all_materials()

        verts_outline = []
        dz_outline = []
        
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
                    created_lenses.append(obj_name)
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
                    elif using_luxcore:
                        materials_bulk, materials_interface = glass_from_Element_luxcore(ele, self.wl, self.mat_refract_only)
                        material_edge = add_blackoutmaterial_luxcore()
                        material_dface = add_diffusematerial_luxcore(viewportcolor=[1, 0, 0, 1])
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
                            if using_cycles or using_luxcore:
                                if i == 0:
                                    paramdict[f'material_name{i+1}'] = materials_bulk[i]
                                elif i == n_loop-1:
                                    paramdict[f'material_name{i+1}'] = materials_bulk[i-1]
                                else:
                                    paramdict[f'material_name{i+1}'] = materials_interface[i-1]
                        for i in range(n_loop-1):
                            paramdict[f'centerthickness{i+1}'] = ele.data['CT'][i0+i]
                        
                        if using_cycles or using_luxcore:
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
                        vo = add_lens(self, context, paramdict=paramdict)
                        verts_outline.append(vo[0]) # only first and last surface
                        verts_outline.append(vo[-1])
                        dz_outline.append(-dz_this)
                        dz_outline.append(-dz_this)
                        #verts_outline = verts_outline + vo
                        #for _ in range(len(vo)):
                        #    dz_outline.append(dz_this)
                        bpy.ops.transform.translate(value=(-dz_this, 0, 0))
                        obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                        created_lenses.append(obj_name)

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
        lens.detector['distance'] = lens.detector['distance'] + self.sensor_offset
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
                obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                created_other.append(obj_name)
            
        t4 = time.perf_counter()
        t41_sum = 0
        t42_sum = 0
        t43_sum = 0

        # add sensor
        sensor_sizex = lens.detector['sizex']/2#*self.sensorfactor
        sensor_sizey = lens.detector['sizey']/2#*self.sensorfactor
        sensorthickness = max(sensor_sizex, sensor_sizey)/20
        zsensor = lens.data['CT_sum'][-1]# + lens.detector['distance']
        if self.addsensor:
            add_sensor(self, context, sensor_sizex, sensor_sizey, zsensor, thicksensor=self.thicksensor, sensorthickness=sensorthickness)
            bpy.ops.transform.translate(value=(-zsensor, 0, 0))
            obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
            created_other.append(obj_name)

        # add camera
        if self.addsensor and self.addcamera:
            # Place camera between sensor and last lens at 10% the spacing.
            # Check first if the sensor distance is unreasonably small, in that case issue a warning
            if lens.detector['distance'] < 0.01*lens.data['CT_sum'][-2]:
                print("[OC] Warning: Distance between last lens and sensor is unreasonably small! The automatically created camera may not yield expected results!")
            zcamera = lens.data['CT_sum'][-2] + 0.9*lens.detector['distance']
            bpy.ops.object.camera_add(location=[-zcamera, 0, 0], rotation=[90*np.pi/180, np.pi, 90*np.pi/180])
            cam = bpy.context.selected_objects[0]
            cam.name = 'OC_Camera'
            cam.data.type = 'ORTHO'
            cam.data.ortho_scale = sensor_sizex*2
            if using_luxcore:
                cam.data.luxcore.imagepipeline.tonemapper.use_autolinear = True
            obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
            created_other.append(obj_name)

        # create housing
        # only makes sense if lenses were created, otherwise verts_outline are empty
        if self.addhousing and self.addlenses:
            if self.housingtype:
                outlinetype='tight'
            else:
                outlinetype='max'
            add_lenshousing_simple(self, context, lens, verts_outline, dz_outline, dshape=self.dshape,
                                   thicksensor=self.thicksensor, sensorthickness=sensorthickness, outlinetype=outlinetype)
            obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
            created_other.append(obj_name)

        # create ray fans    
        if self.addrayfan:
            t40 = time.perf_counter()
            # try to convert the angles and clamp between -90 and 90 degrees
            fanangles = [self.fanangle1]
            try:
                if not self.fanangle_additional in ["", "None", "none", "NONE"]:
                    fanangles_add = [max(min(float(a), 89.99), -89.99) for a in self.fanangle_additional.split(";")]
                    fanangles_add = [angle*np.pi/180 for angle in fanangles_add]
                    fanangles = fanangles_add + fanangles # the "main" fan angle is the last of the list for historic reasons
            except:
                print("FANANGLE ERROR")

            for fanangle in fanangles:
                # set up the rays
                initparams = [self.nrays, self.fandiam*lens.data['rCA'][1], -1*self.fandist, fanangle, self.fanangle2, self.fanangle3] 
                rays = rayfan.RayFan(self.fantype, initparams, store_history=True)
                # try to parse ghost_oder
                ghost_order = self.ghost_order
                nos = lens.num_optical_surfaces
                custom_surflist_valid = False
                if not ghost_order == '':
                    try:
                        surflist = [int(_) for _ in ghost_order.split(',')]
                        if np.any(np.array(surflist) > nos+1) or np.any(np.array(surflist) < 1): 
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
                ocr = self.onlycompleterays

                n_history = len(rays.O_history)
                n_rays = rays.O_history[0].shape[0]
                # create single numpy.array from history. Shape is (n_surf, n_rays, 3)
                O_hist = np.array([rays.O_history[i] for i in range(n_history)])
                idx_complete = ~np.isnan(O_hist[:,:,0].sum(axis=0))

                # reduce array down to fully valid rows
                if ocr:
                    O_hist = O_hist[:, idx_complete, :]
                    n_rays = O_hist.shape[1] # update this variable
                # change axis convention from raytrace to Blender orientation
                O_hist[:, :, [0, 1, 2]] = O_hist[:, :, [2, 0, 1]]
                O_hist[:, :, 0] *= -1
                # always add the verts for the first surface (i.e. ray fan input)
                o1 = np.copy(O_hist[0, :, :]) # initialize first surface
                verts = [list(o1[j, :]) for j in range(n_rays)]
                vidx1 = [j for j in range(n_rays)] # at each surface, the index in the verts list that each point has been assinged
                # loop over pairs of surfaces
                for i in range(1, n_history):
                    offset1 = len(verts)
                    # get new surface
                    o2 = np.copy(O_hist[i, :, :])

                    valid_both = ~np.isnan(o1[:, 0]) & ~np.isnan(o2[:, 0])
                    n_valid = sum(valid_both)
                    valid_any = np.any(valid_both)

                    #add verts and lines
                    if valid_any:
                        verts_new = [list(o2[j, :]) for j, j_valid in enumerate(valid_both) if j_valid]
                        verts = verts + verts_new
                        vidx2 = np.full(n_rays, -1) # invalid points == -1
                        vidx2[valid_both] = range(offset1, offset1 + n_valid)
                        edges_new = [[vidx1[j], vidx2[j]] for j, j_valid in enumerate(valid_both) if j_valid]
                        edges = edges + edges_new
                    else:
                        # can't continue tracing no valid ray
                        break

                    # switch for next loop iteration
                    o2 = np.copy(o1)
                    vidx1 = vidx2
                
                t42 = time.perf_counter()
                t42_sum = t42_sum + t42 - t41
                mesh = bpy.data.meshes.new(name = 'OC_Rayfan')
                mesh.from_pydata(verts, edges, faces)
                obj = object_data_add(context, mesh, operator=self)
                obj_name = obj.name
                created_rayfans.append(obj_name)
                t43 = time.perf_counter()
                t43_sum = t43_sum + t43 - t42

        # create lights
        if self.addlasers:
            created_lights = add_laser_array(self.HFOV, self.VFOV, self.NHFOV, self.NVFOV, self.fandist, lasersize=1.5*2*lens.data['rCA'][1])

        if False: # still the easiest way to block-comment...
            print()
            print('Execution Times:')
            print(f'parsing .zmx-file: {(t1-t0):.3f}')
            print(f'Create Element-Meshes: {(t2-t1):.3f}')
            print(f'Build lens: {(t3-t2):.3f}')
            print(f'Create Aperture-Meshes: {(t4-t3):.3f}')
            if self.addrayfan:
                print(f'Raytracing: {(t41_sum):.3f}')
                print(f'looping ray-fan-mesh: {(t42_sum):.3f}')
                print(f'Adding ray-fan-mesh: {(t43_sum):.3f}')
            print()
        

        # move all components so that the sensor is at the origin
        if self.origin_position == 'sensor':
            bpy.ops.object.select_all(action='DESELECT')
            for objname in created_lenses:
                bpy.data.objects[objname].select_set(True)
            for objname in created_rayfans:
                bpy.data.objects[objname].select_set(True)
            for objname in created_other:
                bpy.data.objects[objname].select_set(True)
            for objname in created_lights:
                bpy.data.objects[objname].select_set(True)
            bpy.ops.transform.translate(value=(zsensor, 0, 0))
        
        # Select a few objects at the end
        # If ray fans were added, select the last one
        # Else, select all lenses
        bpy.ops.object.select_all(action='DESELECT')
        if len(created_rayfans) > 0:
            objname = created_rayfans[-1]
            bpy.data.objects[objname].select_set(True)
        else:
            for objname in created_lenses:
                bpy.data.objects[objname].select_set(True)
        
        if self.display_edit:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='SELECT')

        if not using_cycles and not using_luxcore:
            print('[OC] .zmx-import finished. For transfer to Blender materials,\nuse the following refractive indices')
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
