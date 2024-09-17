import time

import numpy as np

from .import rayfan
from  .import intersect as rt_intersect

import bpy

import warnings
warnings.filterwarnings('ignore')

    

def exec_trace(lens, rays, surfs, trace_detector=True):
    """
    This function performs a trace thorugh all surfaces

    Parameters
    ----------
    lens: function alias
        alias to the lens definition function
    O: np.ndarray
        Origin coordinates of the initial rays
    D: np.ndarray
        Direction vectors of the initial rays
    surfs : iterable
        List of surface indices in order which they shall be traced
    startsurface: int
        Index of the surface, in front of which the original rays start.
        Used in case the tracing shall start at some intermediate point
    trace_detector: bool
        If True, performs the ray tracing to the sensor as the last step.
        If False, an external final trace e.g. to the Blender scene, may be performed

    Returns
    -------
    P_out: list of np.ndarrays
        List of all intermediate hit points of the rays, not including O
    N_out: list of np.ndarrays
        List of all surface normals at the hit points
    OPD_out: TBD
        Optical path travelled by the rays
    P_fail: list of np.ndarrays
        List of hitpoints of the additional rays traced when 

    """
    # initialize variables
    lastsurface = surfs[0] - 1
    n_surfs = len(surfs)
    ld = lens.data # abbreviation

    #start ray tracing loop over surfaces
    for i, idx_s in enumerate(surfs):
        # current direction of propagation
        direction = idx_s - lastsurface
        # check if the next (not current!) surface is identical to the last surface.
        # In this case, a reflection must happen, else a refraction.
        refract = 1
        if i < (n_surfs-1):
            nextsurface = surfs[i+1]
        else:
            nextsurface = max(surfs)+1
        if nextsurface == lastsurface:
            refract = -1

        # get surface parameters
        r = ld['radius'][idx_s]
        rad = ld['rCA'][idx_s]

        # set center of sphere coordinates
        # TODO: Add decenter and tilt
        C = [0, 0, ld['CT_sum'][idx_s] + r]
        C = np.array(C)

        # flip inside/outside IOR depending on direction
        if direction == 1:
            n1, n2 = ld['n'][lastsurface], ld['n'][idx_s]
        else:
            n1, n2 = ld['n'][lastsurface-1], ld['n'][idx_s-1]

        # get the current rays
        O, D = rays.get_rays()
        nfail = sum(np.isnan(O[:,0]))
        
        # Check if aperture is involved
        # ToDo: a more universal way would be to always trace each aperture and check for positive ray-t,
        # comparing with the ray-t of the following lens-surface;
        # That would allow more universal placement but might slow down a lot.
        # ['idx_surface', 'z_ap', 'shape', 'radius', 'n_blades']
        for j, ap in lens.apertures.items():
            aperture_here = False
            if ap['idx_surface'] is not None:
                if ap['idx_surface'] == idx_s-1 and direction == 1:
                    aperture_here = True
                elif ap['idx_surface'] == idx_s and direction == -1:
                    aperture_here = True
                else:
                    continue
            else:
                pass
            if ap['z_ap']:
                zap = ap['z_ap']
            else:
                zap = ld['CT_sum'][idx_s]
            r_ap = ap['radius']
            n_blades = ap['n_blades']
            P, N, idx_fail = rt_intersect.aperture_intersect(O, D, r_ap, zap, n_blades=n_blades)
            rays.update_special_hits(P, idx_fail)
            O[~idx_fail] = float('nan')

        # trace rays
        # calculate lens intersection points
        if r == 0:
            P , N, idx_fail = rt_intersect.circle_intersect(O, D, C[2], rad)
        else:
            P , N, idx_fail = rt_intersect.lens_intersect(O, D, C, r, rad, direction = direction)
        # calculate new ray directions
        if refract == 1:
            D_new = rt_intersect.refract_ray(D, N, n1, n2)
        else:
            D_new = rt_intersect.reflect_ray(D, N)

        # update the origin and direction arrays
        # rays.update(P, D_new, I, OPD, idx_fail)
        a = rays.update(P, D_new, None, None, idx_fail, N=N)
        
        # update parameters for next loop
        lastsurface = int(1*idx_s)

    # TODO: additional aperture check before the detector

    #trace detector
    if trace_detector:    
        O, D = rays.get_rays()
        P, N, idx_fail = rt_intersect.rectangle_intersect(O, D, lens.detector['Quad'])
        rays.update(P, None, None, None, idx_fail, N=N)

    return rays



def trace_to_scene(context, rays):
    EPSILON = 0.001
    P = []
    anyhit = False
    for i in range(rays.O.shape[0]):
        o = np.array(O[i,:]) # + EPSILON
        d = np.array(D[i,:])
        o = o + d*EPSILON
        o[0] = o[0]*-1
        d[0] = d[0]*-1
        scene = context.scene
        graph = bpy.context.evaluated_depsgraph_get()
        result = scene.ray_cast(graph, origin=o, direction = d)
        hit, p, normal, index, ob, matrix = result
        if hit:
            anyhit = True
            p[0] = p[0]*-1
            P.append(p)
        else:
            P.append([float('nan'), float('nan'), float('nan')])
    P = np.array(P)
    idx_fail = np.isnan(P[:,0])
    rays.update(P, None, None, None, idx_fail)
    
    return rays