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

import time

import numpy as np

from . import intersect as rt_intersect
from . import intersect_asphere as rt_intersect_asph
from ..surface.toric import get_z_toric, get_N_toric

import warnings
warnings.filterwarnings('ignore')


def exec_trace(lens, rays, surfs=None, trace_detector=True):
    """
    This function performs a trace thorugh all surfaces

    Parameters
    ----------
    lens: LensSystem() instance
        The Lens system to be traced
    rays: Rayfan() instance
        A Rayfan-class object containing the initial rays for tracing
    surfs : iterable
        List of surface indices in order which they shall be traced
    trace_detector: bool
        If True, performs the ray tracing to the sensor as the last step.
        If False, allows an external final trace, e.g. to the Blender scene

    Returns
    -------
    rays: Rayfan() instance
        The updated Rayfan object containing the final ray positions

    """
    # initialize variables
    if lens.surf_sequence is not None:
        surfs = lens.surf_sequence
    elif surfs is None:
        # default to just tracing straight through the surfaces
        surfs = [i+1 for i in range(lens.num_surfaces)]
    lastsurface = surfs[0] - 1
    n_surfs = len(surfs)
    ld = lens.data # abbreviation
    reflectionstate = 1 # TODO: merge this into a better direction handling mechanism

    # start ray tracing loop over surfaces
    for i, idx_s in enumerate(surfs):
        # current direction of propagation
        direction = int(np.sign(idx_s - lastsurface))
        # check if the next (not current!) surface is identical to the last surface.
        # In this case, a reflection must happen, else a refraction.
        refract = 1
        ismirror = ld['ismirror'][idx_s]
        ismirror_pre = ld['ismirror'][idx_s-1]
        mirrorfactor = 1 - 2*ismirror_pre
        reflectionstate = reflectionstate*mirrorfactor
        if i < (n_surfs-1):
            nextsurface = surfs[i+1]
        else:
            nextsurface = max(surfs)+1
        if nextsurface == lastsurface:
            # case ghost
            refract = -1
        if ismirror:
            # case mirror
            refract = -1

        # get surface parameters
        rad = ld['rCA'][idx_s]
        r = ld['radius'][idx_s]
        k = ld['asph'][idx_s][0]
        A = ld['asph'][idx_s][1:]
        r2 = ld['radius2'][idx_s]
        k2 = ld['asph2'][idx_s][0]
        A2 = ld['asph2'][idx_s][1:]
        # determine the type of surface w.r.t. the intersection algorithms:
        surftype = ld['type'][idx_s]
        # surface rotation for e.g. cylinder lenses
        surf_rotation = ld['surf_rotation'][idx_s]

        # set center of sphere coordinates
        # TODO: Add decenter and tilt
        C = [0, 0, ld['CT_sum'][idx_s] + r]
        C = np.array(C)
        C_CT = [0, 0, ld['CT_sum'][idx_s]]
        C_CT = np.array(C_CT)

        # flip inside/outside IOR depending on direction
        if direction == 1:
            n1, n2 = ld['n'][lastsurface], ld['n'][idx_s]
        else:
            n1, n2 = ld['n'][lastsurface-1], ld['n'][idx_s-1]

        # get the current rays
        O, D = rays.get_rays()
        
        # Check if aperture is involved
        # A more universal way would be to always trace each aperture and check for positive ray-t,
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
                # TODO: Check if aperture is between previous and next surface to save unneccessary tracing
                aperture_here=True
            if ap['z_ap']:
                zap = ap['z_ap']
            else:
                zap = ld['CT_sum'][idx_s]
            if aperture_here:
                r_ap = ap['radius']
                n_blades = ap['n_blades']
                P_inside, P_outside, N, idx_fail_i, idx_fail_o = rt_intersect.aperture_intersect(O, D, r_ap, zap, n_blades=n_blades)
                rays.update_special_hits(P_outside, idx_fail_i)
                O[~idx_fail_o] = float('nan')

        # calculate lens intersection points
        if surftype == 'aspheric':
            P, N, idx_fail = rt_intersect_asph.intersect_asphere(O, D, C_CT, rad, r, k, A)
        elif surftype == 'toric':
            z_fun_params = [r, r2, surf_rotation]
            P, N, idx_fail = rt_intersect_asph.intersect_implicit(O, D, C_CT, rad, get_z_toric, z_fun_params,
                                                                  N_fun=get_N_toric)
        else:
            P, N, idx_fail = rt_intersect.lens_intersect(O, D, C_CT, r, rad,
            k=k, A=A, surf_rotation=surf_rotation, surfshape=surftype, direction=direction)#*reflectionstate)
        
        # calculate new ray directions
        if refract == 1:
            D_new = rt_intersect.refract_ray(D, N, n1, n2)
        else:
            D_new = rt_intersect.reflect_ray(D, N)

        # update the origin and direction arrays
        rays.update(P, D_new, None, None, idx_fail, N=N)
        
        # update parameters for next loop
        lastsurface = int(1*idx_s)

        if np.all(idx_fail):
            print('ALL FAIL at surface', idx_s)
            return rays

    # additional aperture check before the detector
    O, D = rays.get_rays()
    for j, ap in lens.apertures.items():
        aperture_here = False
        if ap['idx_surface'] is not None:
            if idx_s == ap['idx_surface']: # idx_s has retained the value of the last surface from the loop
                aperture_here = True
            else:
                continue
        else:
            # TODO: Check if aperture is between previous and next surface to save unneccessary tracing
            aperture_here=True
        if ap['z_ap']:
            zap = ap['z_ap']
        else:
            zap = ld['CT_sum'][idx_s]
        if aperture_here:
            r_ap = ap['radius']
            n_blades = ap['n_blades']
            P_inside, P_outside, N, idx_fail_i, idx_fail_o  = rt_intersect.aperture_intersect(O, D, r_ap, zap, n_blades=n_blades, pass_inside=True)
            rays.update_special_hits(P_outside, idx_fail_o)
            # reset the origin rays
            rays.update(O, D, None, None, idx_fail_i, N=N) # use idx_fail_i because aperture_intersect has opposite logic (aperture IS hit)

    # trace detector
    if trace_detector:
        O, D = rays.get_rays()
        P, N, idx_fail = rt_intersect.rectangle_intersect(O, D, lens.detector['Quad'])
        rays.D_tosensor = np.array(rays.D)
        rays.update(P, None, None, None, idx_fail, N=N)

    return rays