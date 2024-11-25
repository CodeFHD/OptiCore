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

import re
import os
import time
import numpy as np

from mathutils import Vector

from ..utils.raytrace.element import Element
from ..utils.raytrace.lenssystem import Lenssystem
from ..utils.raytrace.trace_sequential import exec_trace, trace_to_scene
from ..utils.raytrace import rayfan
from ..elements import add_lens, add_mirror, get_default_paramdict_lens, get_default_paramdict_mirror, add_circular_aperture

import bpy
from bpy.props import FloatProperty, IntProperty, EnumProperty, StringProperty, BoolProperty, FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add


DEBUGPRINT = False
def debugprint(*messages):
    if len(messages) == 0: messages = ''
    if DEBUGPRINT: print(*messages)


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
                    ("3D_square","3D square",""),
                    ("3D_tri","3D tris",""),
                    ("3D_random","3D random",""),
                    ("2D_finite","2D Finite",""),
                    ("3D_rings", "3D rings", ""),},
           default = "2D",
           description="Ray Fan Type",
           #options={'HIDDEN'},
           )
    fandist : FloatProperty(
           name="Fan Origin",
           default = 20.,
           description="Distance where ray fan originates.",
           )
    fanangle : FloatProperty(
           name="Fan Angle",
           default = 0.,
           description="Angle of Ray Fan.",
           min = -90.,
           max = 90.,
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
        col.prop(self, 'fanangle')
        col.prop(self, 'fanangle_additional')
        col.prop(self, 'fandiam')
        col.prop(self, 'aperturediam')
        col.prop(self, 'excludedetector')
        col.prop(self, 'addlenses')
        col.prop(self, 'addaperture')
        col.prop(self, 'addrayfan')
        col.prop(self, 'tracetoscene')
        col.prop(self, 'autofocus')
        col.prop(self, 'onlycompleterays')
        col.label(text="Mechanical Parameters")
        col.prop(self, 'dshape')
        col.prop(self, 'display_edit')
        col.prop(self, 'split_cemented')
        col.prop(self, 'num1')
        col.prop(self, 'num2')

    
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
        lens = load_from_zmx(self.path)
        lens.apertures[0]['radius'] = lens.apertures[0]['radius']*self.aperturediam
        created_objects = [] # list of the individual objects that were created
        t1 = time.perf_counter()
        
        # create meshes from the elements
        if self.addlenses:
            CT_sum = 0
            radius_outer = 0
            for ele, dz in lens.elements:
                """
                TODO: Possibly to be foreseen in this implementation:
                In case of a multi-bounce system, i.e. where the same surface is used multiple times,
                implement a logic that it is not necessary to create the same element twice.
                """
                num_surfaces = len(ele.data['radius'])
                if num_surfaces > 4 and not self.split_cemented:
                    # TODO: Currently only supporting up to triplet groups
                    continue
                dz_this = CT_sum + dz
                CT_sum = CT_sum + dz + sum(ele.data['CT'][:-1])
                if ele.ismirror:
                    paramdict = get_default_paramdict_mirror()
                    paramdict['rad'] = -1*ele.data['radius'][0]*ele.direction
                    paramdict['num1'] = self.num1
                    paramdict['num2'] = self.num2
                    paramdict['mirrorradius'] = ele.data['rCA'][0]
                    paramdict['centerthickness'] = 1.
                    paramdict['mtype'] = 'aspheric'
                    paramdict['k'] = ele.data['asph'][0][0]
                    paramdict['A'] = ele.data['asph'][0][1:]
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
                elif self.split_cemented and num_surfaces > 2:
                    for i in range(num_surfaces - 1):
                        paramdict = get_default_paramdict_lens()
                        paramdict[f'rad1'] = ele.data['radius'][i]
                        paramdict[f'rad2'] = ele.data['radius'][i+1]
                        paramdict['num1'] = self.num1
                        paramdict['num2'] = self.num2
                        paramdict['lensradius'] = max(ele.data['rCA'][i:i+2])
                        if ele.data['rCA'][i] < max(ele.data['rCA'][i:i+2]):
                            paramdict[f'flangerad1'] = max(ele.data['rCA'][i:i+2]) - ele.data['rCA'][i]
                        if ele.data['rCA'][i+1] < max(ele.data['rCA'][i:i+2]):
                            paramdict[f'flangerad2'] = max(ele.data['rCA'][i:i+2]) - ele.data['rCA'][i+1]
                        paramdict[f'centerthickness1'] = ele.data['CT'][i]
                        test_asph = [x in [0, None] for x in ele.data['asph'][i]]
                        if not np.all(test_asph):
                            paramdict[f'k1'] = ele.data['asph'][i][0]
                            paramdict[f'A1'] = ele.data['asph'][i][1:]
                            paramdict[f'ltype1'] = 'aspheric'
                        test_asph = [x in [0, None] for x in ele.data['asph'][i+1]]
                        if not np.all(test_asph):
                            paramdict[f'k2'] = ele.data['asph'][i+1][0]
                            paramdict[f'A2'] = ele.data['asph'][i+1][1:]
                            paramdict[f'ltype2'] = 'aspheric'
                        paramdict['makedoublet'] = '1'
                        paramdict['dshape'] = self.dshape
                        add_lens(self, context, paramdict=paramdict)
                        bpy.ops.transform.translate(value=(-dz_this, 0, 0))
                        obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                        created_objects.append(obj_name)
                        dz_this = dz_this + ele.data['CT'][i]
                else:
                    paramdict = get_default_paramdict_lens()
                    for i in range(num_surfaces):
                        paramdict[f'rad{i+1}'] = ele.data['radius'][i]
                    paramdict['num1'] = self.num1
                    paramdict['num2'] = self.num2
                    paramdict['lensradius'] = max(ele.data['rCA'])
                    for i in range(num_surfaces):
                        if ele.data['rCA'][i] < max(ele.data['rCA']):
                            paramdict[f'flangerad{i+1}'] = max(ele.data['rCA']) - ele.data['rCA'][i]
                    for i in range(num_surfaces-1):
                        paramdict[f'centerthickness{i+1}'] = ele.data['CT'][i]
                    for i in range(num_surfaces):
                        test_asph = [x in [0, None] for x in ele.data['asph'][i]]
                        if not np.all(test_asph):
                            paramdict[f'k{i+1}'] = ele.data['asph'][i][0]
                            paramdict[f'A{i+1}'] = ele.data['asph'][i][1:]
                            paramdict[f'ltype{i+1}'] = 'aspheric'
                    #if num_surfaces == 3:    
                    paramdict['makedoublet'] = str(num_surfaces - 1)
                    paramdict['dshape'] = self.dshape
                    # paramdict['optiverts'] = 
                    # paramdict['material_name'] = 
                    # paramdict['material_name2'] = 
                    # paramdict['material_name3'] = 
                    # paramdict['shade_smooth'] = 
                    # paramdict['smooth_type'] = 
                    # paramdict['ior'] = 
                    # paramdict['display_edit'] =
                    add_lens(self, context, paramdict=paramdict)
                    bpy.ops.transform.translate(value=(-dz_this, 0, 0))
                    obj_name = bpy.context.selected_objects[0].name # assuming only one is selected
                    created_objects.append(obj_name)
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
            fanangles = [self.fanangle]
            try:
                if not self.fanangle_additional in ["", "None", "none", "NONE"]:
                    fanangles_add = [max(min(float(a), 89.99), -89.99) for a in self.fanangle_additional.split(";")]
                    fanangles = fanangles_add + fanangles
            except:
                print("FANANGLE ERROR")

            for fanangle in fanangles:
                # set up the rays
                # TODO: it might be a less branched approach to set up the initparams with a dictionary. Consider with future updates
                if self.fantype in ["2D_finite"]:
                    rayfany = self.fandist*np.tan(fanangle*np.pi/180)
                    initparams = [self.nrays, lens.data['rCA'][1]*self.fandiam, -1*self.fandist, -1.*rayfany]
                else:# self.fantype in ["2D"]:
                    initparams = [self.nrays, lens.data['rCA'][1]*self.fandiam, -1*self.fandist, fanangle*np.pi/180] 
                #else:
                #    initparams = [self.nrays, lens.data['rCA'][1]*self.fandiam, -1*self.fandist, self.fanangle*np.pi/180]    
                rays = rayfan.RayFan(self.fantype, initparams, store_history=True)
                surflist = [i for i in range(1, lens.num_optical_surfaces + 1)] # standard surface list for non-ghost trace
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
        
        # Select all just created objects
        bpy.ops.object.select_all(action='DESELECT')
        for o in bpy.data.objects:
            if o.name in created_objects:
                o.select_set(True)
        
        if self.display_edit:
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action='SELECT')
        
        return {'FINISHED'}

def split_file(f):
    lines = f.readlines()
    headerlines = []
    surflines = {} # dict, for each surface index, a list of the lines
    footerlines = []
    lc = 0 # linecounter

    # loop over header
    for line in lines:
        line = line.strip()
        if line.startswith('SURF'):
            break
        headerlines.append(line)
        lc += 1

    # loop over surfaces
    surfidx = -1
    for line in lines[lc:]:
        if re.match(r'\w', line) and not line.startswith('SURF'):
            # break when something else than a surf comes up
            break
        line = line.strip()
        lc += 1
        if line.startswith('SURF'):
            if surfidx > -1:
                surflines[surfidx] = curlines 
            surfidx = int(line.split()[1])
            curlines = []
        else:
            curlines.append(line)
    # finish last surface
    surflines[surfidx] = curlines 

    # loop over footer
    for line in lines[lc:]:
        line = line.strip()
        footerlines.append(line)

    return headerlines, surflines, footerlines


def parse_surface(surflines):
    surf_info = {} # store in one dictionary for not passing too many variables
    # defaults if parameters not given for this surface
    type = None
    isstop = False
    curv = None
    CT = None
    hasglass = False
    ismirror = False
    glass = None
    rCA = None
    lrad = None
    k = 0
    A = []
    # search the parameters
    for line in surflines:
        if line.startswith('TYPE'):
            type = line.split()[1]
            if type == 'XASPHERE': type = 'EVENASPH'
        elif line.startswith('STOP'):
            isstop = True
        elif line.startswith('CURV'):
            curv = float(line.split()[1])
        elif line.startswith('DISZ'):
            CT = float(line.split()[1])
        elif line.startswith('GLAS'):
            hasglass = True
            glass = line.split()[1]
            if glass == '___BLANK':
                glass = ' '.join(line.split()[1:])
            if glass == 'MIRROR':
                ismirror = True
                hasglass = False
        elif line.startswith('DIAM'):
            rCA = float(line.split()[1])
        elif line.startswith('MEMA'):
            lrad = float(line.split()[1])
        elif line.startswith('CONI'):
            k = float(line.split()[1])
        elif line.startswith('PARM'):
            iparm = int(line.split()[1])
            if type == 'EVENASPH':
                if iparm == 1:
                    # I have yet to see a file that has a power-2 coefficient
                    # (assuming iparm == 1 is for power 2 coefficient)
                    # Also the geometry function currently starts at power-4
                    pass
                else:
                    A.append(float(line.split()[2]))
        # TODO: parse aspheric surface
    # fill dict and return
    if A == []: A = [0]
    surf_info['type'] = type
    surf_info['isstop'] = isstop
    surf_info['radius'] = 0 if curv == 0 else 1./curv
    surf_info['CT'] = CT
    surf_info['hasglass'] = hasglass
    surf_info['ismirror'] = ismirror
    surf_info['glass'] = glass
    surf_info['rCA'] = rCA
    surf_info['lrad'] = lrad
    surf_info['asph'] = [k] + A
    return surf_info

# Because of FLoat rounding, values are rounded to these digits for comparison
EPSILON_DIGITS = 3

def identify_elements(surf_infos):
    """
    This function evaluates the previously parsed surfaces and splits it into physical Elements,
    i.e. cemented pieces of glass, or mirrors.
    - Dummy surfaces are discarded (except their CT is counted).
    - Double passes thorugh one Element are detected.
    - The order in which the so created filtered surfaces are passed is determined.
    """
    idx_elements = [] # one list for each element, with the corresponding surface indices
    idx_order = [] # lists the order that surfaces were passed, in case of reflections and multiple passes
    cur_element = []
    direction = True # direction is flipped at mirrors
    
    # some varaible for tracking duplicate surfaces due to reflections, multiple passes, or split-surface-tables
    abspos = 0 # position of the current surface in absolute units
    cur_idx = 0 # count of unique surfaces
    abspos_at_idx = [] # for the unique surfaces: position along the optical axis
    radius_at_idx = [] # for the unique surfaces: radius of curvature
    glass_at_idx = [] # for the unique surfaces: if there is a glass or mirror at that index
        
    # for detection of multi-surface splits, keep track of the last proper surface radius
    abspos_last = -999
    absrad_last = 0
            
    CT_cumulative = {} # cumulative center thickness including dummies. 
    idx_lastvalid = -1
    
    for i in range(1, len(surf_infos)-1): # loop this way just to make sure surfaces are ordered and skip surfaces 0 and -1 (object and image)
        
        debugprint()
        debugprint('IDENTIFY', i)
        debugprint('-'*40)

        glass_prev = surf_infos[i-1]['hasglass']
        glass_here = surf_infos[i+0]['hasglass']
        glass_next = surf_infos[i+1]['hasglass']
        mirror_prev = surf_infos[i-1]['ismirror']
        mirror_here = surf_infos[i+0]['ismirror']
        mirror_next = surf_infos[i+1]['ismirror']
        glasstype_prev = surf_infos[i-1]['glass']
        glasstype_here = surf_infos[i+0]['glass']
        glasstype_next = surf_infos[i+1]['glass']
        
        debugprint('GLASSES', glass_prev, glass_here, glass_next)
        debugprint('GLASSTYPES', glasstype_prev, glasstype_here, glasstype_next)
        debugprint('MIRRORS', mirror_prev, mirror_here, mirror_next)
        
        absrad_prev = round(abs(surf_infos[i-1]['radius']), EPSILON_DIGITS)
        absrad_here = round(abs(surf_infos[i+0]['radius']), EPSILON_DIGITS)
        absrad_next = round(abs(surf_infos[i+1]['radius']), EPSILON_DIGITS)
        
        abspos_here = round(abspos, EPSILON_DIGITS) # abspos needs to be incremented whatever happens, but the value before increment counts
        abspos = abspos + surf_infos[i]['CT'] # carry without rounding
        abspos_r = round(abspos, EPSILON_DIGITS)
                          
        isdummy = ((glasstype_here == glasstype_prev) or (not glasstype_here and not glass_prev and cur_element == [])) and not mirror_here
        # case 1: The surface is a dummy
        # no glass, nor refraction; no mirror, no reflection
        if isdummy:
            debugprint('DUMMY')
            debugprint(abspos_here, absrad_here, abspos_last, absrad_last, glass_here)
            abspos_last = abspos_here
            if idx_lastvalid != -1:
                CT_cumulative[idx_lastvalid] = CT_cumulative[idx_lastvalid] + surf_infos[i]['CT'] 
            # nothing to do for dummy surfaces
            continue
        # when the dummy check is passed, create the new CT entry
        idx_lastvalid = i
        CT_cumulative[i] = surf_infos[i]['CT'] 
        
        # isopticalsurf = glass_here or glass_prev or mirror_here # condition that the surface is optical
        lastopticalsurf = (glass_prev and not glass_here and not mirror_prev) or mirror_here # condition that the surface ends a stack of optical surfaces. TODO: This assumes that mirrors always end the element, half-mirrors (e.g. polarized pancake lenses) could break this.
        # Check if the following surface would be a duplicate to this one
        split_next = (abspos_r == abspos_here) and (absrad_next == absrad_here)# and glass_here
        split_continued = (abspos_here == abspos_last) and (absrad_here == absrad_last) and i>1
        # Duplicate surface is defined by vertex position along optical axis and absolute radius of curvature
        # TODO, could keep track of reflections and consider sign of radius
        duplicate_here = (abspos_here in abspos_at_idx) and (absrad_here in radius_at_idx) and not split_continued
        
        debugprint(abspos_here, absrad_here, abspos_last, absrad_last, glass_here)
        
        # case 2: The surface is a duplicate of a previous surface.
        if duplicate_here:
            prev_idx = abspos_at_idx.index(abspos_here)
            debugprint('DUPLICATE', prev_idx)
            # case 2.1: The duplicate is replaced because it features glass
            if not glass_at_idx[prev_idx] and direction*glass_here:
                cur_element[-1] = cur_element[-1] + 1 # this condition should only be invoked at the first pass of this surface pair so this should be safe
                glass_at_idx[prev_idx] = True
            # case 2.2: The previous surface is taken as-is
            else:
                idx_order.append(prev_idx + 1) # the tracer later counts starting at surface 1
            abspos_last = abspos_at_idx[prev_idx]
            absrad_last = radius_at_idx[prev_idx]
        # case 3: The surface is new
        else:
            debugprint('NEW', cur_idx)
            # case 3.1: check if the surface is a continuation of a split. If not, finish previous element and start new
            if split_continued:
                debugprint('CONTINUING SPLIT')
                cur_element[-1] = i
                glass_at_idx[-1] = glass_here
                continue
            cur_element.append(i)
            abspos_at_idx.append(abspos_here)
            radius_at_idx.append(absrad_here)
            glass_at_idx.append(glass_here or mirror_here)
            idx_order.append(cur_idx + 1) # the tracer later counts starting at surface 1
            cur_idx = cur_idx + 1
            abspos_last = abspos_here
            absrad_last = absrad_here
            debugprint('LEN', len(abspos_at_idx))
            # case 3.1: It finishes an Element
            if lastopticalsurf and not split_next:
                debugprint('FINISHING', cur_element)
                idx_elements.append(cur_element)
                cur_element = []    
            # case 3.2: It does not finish an Element
            #else:
            #    pass

        if mirror_here:
            direction = not direction
           
    # If there is an unfinished element for the last surface, termiante it
    if cur_element != []:
        debugprint('FINISHING AT END', cur_element)
        idx_elements.append(cur_element)
        
    return idx_elements, idx_order, CT_cumulative


def get_stop(surf_infos, idx_first=0):
    stopidx = 1
    for idx, surf in surf_infos.items():
        if idx == 1:
            if surf['rCA'] is not None:
                stoprad = surf['rCA']
            elif surf['lrad'] is not None:
                stoprad = surf['lrad']
            break
    for idx, surf in surf_infos.items():
        if surf['isstop']:
            stopidx = idx
            if surf['rCA'] is not None:
                stoprad = surf['rCA']
            elif surf['lrad'] is not None:
                stoprad = surf['lrad']
    if idx_first > stopidx:
        CT_list = [surf_infos[i]['CT'] for i in range(stopidx, idx_first)]
        z_stop = -sum(CT_list)
    else:
        CT_list = [surf_infos[i]['CT'] for i in range(idx_first, stopidx)]
        z_stop = sum(CT_list)
    return stopidx, stoprad, z_stop

def create_elements(surf_infos, idx_elements, CT_cumulative):
    elements = []
    # firstElement = True
    direction = 1 # at a mirror, the direction is flipped and CT is negative
    for i, idx_list in enumerate(idx_elements):
        ele = Element()
        if len(idx_list) == 1:
            if surf_infos[idx_list[0]]['glass'] == 'MIRROR':
                ele.ismirror = True
        # Fill the element data
        for j, idx in enumerate(idx_list):
            ismirror = surf_infos[idx]['glass'] == 'MIRROR'
            ele.data['type'].append(surf_infos[idx]['type'])
            ele.data['radius'].append(surf_infos[idx]['radius'])
            ele.data['asph'].append(surf_infos[idx]['asph'])
            # ele.data['lrad'].append(surf_infos[idx]['lrad'])
            ele.data['rCA'].append(surf_infos[idx]['rCA'])
            ele.data['surf_decenter'].append([0,0]) # (surf_infos[idx][])
            ele.data['surf_tilt'].append([0,0]) # (surf_infos[idx][])
            if ismirror:
                ele.data['coating'].append('MIRROR')
            else:
                ele.data['coating'].append(None)
            if True:#j < len(idx_list) - 1 or ele.ismirror:
                ele.data['CT'].append(CT_cumulative[idx]*direction)# (surf_infos[idx]['CT']*direction)
                ele.data['material'].append(surf_infos[idx]['glass'])
        # get the element offsets
        if i == 0:
            dz = 0
        else:
            CT_list = [surf_infos[j]['CT'] for j in range(idx_elements[i-1][-1], idx_list[0])]
            dz = sum(CT_list)
        #if direction == -1:
        #    ele.reverse()
        ele.direction = direction    
        if ele.ismirror:
            direction = -1*direction
        elements.append([ele, dz])
        debugprint()
        debugprint('FINISHED ELEMENT')
        debugprint()
    # get the sensor distance
    CT_list = [surf_infos[j]['CT'] for j in range(idx_elements[-1][-1], len(surf_infos))]
    dz_sensor = sum(CT_list)
    return elements, dz_sensor

def load_from_zmx(filename, testmode=False, logfile=None):
    # check encoding of file
    f = open(filename, 'rb')
    data = f.read(2)
    if data == b'\xff\xfe':
        encoding = 'utf-16-le'
    else:
        encoding = 'utf-8'
    f.close()
    #with open (filename, 'r', encoding='utf-16-le') as f:
    with open (filename, 'r', encoding=encoding) as f:    
        he, su, fo = split_file(f)
        surf_infos = {}
    for idx, surface in su.items():
        surf_info = parse_surface(surface)
        surf_infos[idx] = surf_info
    idx_elements, idx_order, CT_cumulative = identify_elements(surf_infos)
    debugprint()
    debugprint('idx order', idx_order)
    debugprint()
    elements, dz_sensor = create_elements(surf_infos, idx_elements, CT_cumulative)

    lens = Lenssystem()
    lens.testmode = testmode
    lens.logfile = logfile
    lens.elements = elements
    lens.surf_sequence = idx_order

    lens.num_optical_surfaces = sum([len(l) for l in idx_elements])
    
    idx_imagesurf = len(surf_infos)-1
    if surf_infos[idx_imagesurf]['hasglass']:
        lens.imagemat = surf_infos[idx_imagesurf]['glass']
    
    pixelpitch = 0.01 # default because that is not stored in zmx files AFAIK
    if surf_infos[idx_imagesurf]['rCA'] is not None:
        imsize = int(np.ceil(surf_infos[idx_imagesurf]['rCA']/pixelpitch))
        # multiply by 10 because usually the zmx file will specify the image size up to a design angle
        # but not including vignetting. And full vignetting can exceed it by far.
        # This doesn't affect performance, only ray fan cutoff
        lens.detector['sizex'] = imsize*10
        lens.detector['sizey'] = imsize*10
    else:
        lens.detector['sizex'] = 10240 # default
        lens.detector['sizey'] = 10240 # default
    lens.detector['distance'] = dz_sensor
    lens.detector['pixelpitch'] = pixelpitch # default

    stopidx, stoprad, z_stop = get_stop(surf_infos, idx_elements[0][0])
    lens.apertures[0] = {'idx_surface': None, 'z_ap': z_stop, 'shape': 'CIRCLE', 'radius': stoprad, 'n_blades': 0}

    if testmode:
        lens.detector['distance'] = 5000
        # build the lens, raytrace, and get the rays for plotting
        lens.build(0.586)

    return lens