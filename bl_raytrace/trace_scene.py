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

def trace_to_scene(context, rays):
    EPSILON = 0.001
    P = []
    anyhit = False
    O, D = rays.get_rays()
    for i in range(O.shape[0]):
        o = np.array(O[i,:])
        d = np.array(D[i,:])
        o = o + d*EPSILON
        o[[0,1,2]] = o[[2,0,1]]
        d[[0,1,2]] = d[[2,0,1]]
        o[0] = o[0]*-1
        d[0] = d[0]*-1
        scene = context.scene
        graph = bpy.context.evaluated_depsgraph_get()
        result = scene.ray_cast(graph, origin=o, direction = d)
        hit, p, normal, index, ob, matrix = result
        p = np.array(p)
        p[[2,0,1]] = p[[0,1,2]]
        if hit:
            anyhit = True
            p[2] = p[2]*-1
            P.append(p)
        else:
            P.append([float('nan'), float('nan'), float('nan')])
    P = np.array(P)
    idx_fail = np.isnan(P[:,0])
    rays.update(P, None, None, None, idx_fail)
    
    return rays