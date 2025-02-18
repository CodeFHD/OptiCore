import os
import numpy as np

import bpy
from bpy_extras.object_utils import AddObjectHelper

from . import load_zmx


"""
THIS FILE IS NOT COMMITTED
FUNCTIONS MUST BE DIABLED IN MASTER __init__.py
BEFORE PUSHING A COMMIT!
"""


basefolder = r'C:\Users\johannes\Documents\0_optics\lensdesigns_com'
design_types = ['photo_prime', 'photo_zoom', 'telescope', 'eyepiece', 'microscope', 'projector']

class OBJECT_OT_test_zmx(bpy.types.Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.test_zmx"
    bl_label = "Test zmx files"
    bl_options = {'REGISTER', 'UNDO'}
    
    def draw(self, context):
        pass
    
    def execute(self, context):
        test_zemax_library()
        postprocess_valid_lenses()
        return {'FINISHED'}

def test_zemax_library():
    for design_type in design_types:
        folder = os.path.join(basefolder, design_type)
        logfile = os.path.join(basefolder, 'log_' + design_type + '.txt')
        f_log = open(logfile, 'a')
        
        design_files = next(os.walk(folder))[2]
        zmx_files = []
        for f in design_files:
            if 'nonseq' in f.lower():
                continue
            if f.lower().endswith('.zmx'): zmx_files.append(f)
        
        for f in zmx_files:
            print()
            print('NOW: ', f)
            fname = os.path.join(folder, f)
            f_log.write(f)
            try:
                lens = load_zmx.load_from_zmx(fname, testmode=True, logfile=f_log)
            except:
                f_log.write(' FAILED')
            f_log.write('\n')

        f_log.close()
        
def postprocess_valid_lenses():
    usablefile = os.path.join(basefolder, 'functional_lenses.txt')
    f_use = open(usablefile, 'a')
    glassfile = os.path.join(basefolder, 'missing_glasses.txt')
    f_glas = open(glassfile, 'a')
    glassfile2 = os.path.join(basefolder, 'missing_glasses_rii.txt')
    f_glas2 = open(glassfile2, 'a')
    riifile = r'C:\Users\johannes\Documents\software\py\py_simpletrace\data\refractiveindexinfo\database\data-nk\rii_mapping.dat'
    f_rii = open(riifile)
    glasses_rii = []
    glasses = []
    for line in f_rii.readlines():
        glasscode = line.split()[0]
        glasses_rii.append(glasscode.lower())
    f_rii.close()
    for design_type in design_types:
        f_use.write(design_type + '\n')
        f_use.write('----------------------------------------------\n')
        logfile = os.path.join(basefolder, 'log_' + design_type + '.txt')
        f_log = open(logfile, 'r')
        for line in f_log.readlines():
            n = len(line.split())
            if 'FAILED' in line:
                continue
            elif n==1:
                f_use.write(line)
            else:
                glasslist = line.split()[1:]
                for g in glasslist:
                    if g == 'FAILED':
                        pass
                    elif not g in glasses:
                        glasses.append(g)
        f_use.write('\n')
        f_log.close()
    f_use.close()
    glasses = sorted(glasses)
    l1 = glasses[0][0]
    for g in glasses:
        l2 = g[0]
        if l2 != l1:
            f_glas.write('\n')
            f_glas2.write('\n')
            l1 = l2
        f_glas.write(g + ' ')
        if g.lower() not in glasses_rii:
            f_glas2.write(g + ' ')
            
        
    f_glas.close()
    f_glas2.close()

if __name__ == '__main__':
    #test_zemax_library()
    postprocess_valid_lenses()

