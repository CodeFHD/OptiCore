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

from .flat import *
from .parabolic import *
from .spherical import *
from .aspheric import *
from .rim import *
from .toric import *
from .sagsurface import *


def add_surface(ltype, surf_subtype, squarelens, N1, N2, lrad,
                srad, k=None, A=None, RY2=None, k2=None, A2=None,
                surf_rotation=0,
                zadd=0, nVerts=0,
                dshape=False, lrad_ext=None):
    """
    This function is a wrapper for the different surface shapes
    and returns the corresponding vertices, surface index lists and normals

    to keep the code short observe the following variable names:
    v, f, n == verts, faces, normals
    """

    """Rectangular cutout"""
    if squarelens:
        if ltype == 'flat' or surf_subtype == 'flat':
            # flat surface
            v, f, n, v_outline = add_sqflat_surface(lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape)
        elif ltype == 'rotational' and surf_subtype == 'spherical':
            # normal spherical surface
            v, f, n, v_outline = add_sqspherical_surface(srad, lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape, lwidth_ext=lrad_ext)
        elif ltype == 'rotational' and surf_subtype != 'spherical':
            # normal aspheric surface
            v, f, n, v_outline = add_sagsurface_rectangular(srad, k, A, lrad, N1, N2, surf_rotation=surf_rotation, zadd=zadd, nVerts=nVerts, surftype='aspheric', dshape=dshape, lrad_ext=lrad_ext)
        elif ltype == 'cylindrical' and surf_subtype == 'spherical':
            # normal aspheric surface
            v, f, n, v_outline = add_sagsurface_rectangular(srad, k, A, lrad, N1, N2, surf_rotation=surf_rotation, zadd=zadd, nVerts=nVerts, surftype='cylindrical', dshape=dshape, lrad_ext=lrad_ext)
        elif ltype == 'cylindrical' and surf_subtype != 'spherical':
            # normal aspheric surface
            v, f, n, v_outline = add_sagsurface_rectangular(srad, k, A, lrad, N1, N2, surf_rotation=surf_rotation, zadd=zadd, nVerts=nVerts, surftype='acylindrical', dshape=dshape, lrad_ext=lrad_ext)
            """
        elif ltype == 'cylindricX' and surf_subtype == 'spherical':
            # X-cylinder spherical
            v, f, n = sfc.add_sqspherical_surface(srad, lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape, lwidth_ext=lrad_ext, cylinderaxis='X')
        elif ltype == 'cylindricY' and surf_subtype == 'spherical':
            # Y-cylinder spherical
            v, f, n = sfc.add_sqspherical_surface(srad, lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape, lwidth_ext=lrad_ext, cylinderaxis='Y')
            """
        elif ltype == 'toric' and surf_subtype == 'spherical':
            v, f, n, v_outline = add_sagsurface_rectangular(srad, k, A, lrad, N1, N2, zadd=zadd, nVerts=nVerts,
                    rad2=RY2, k2=None, A2=None,
                    surf_rotation=surf_rotation, surftype='toric',dshape=dshape, lrad_ext=lrad_ext)
        else:
            # in case of anything not covered by the above
            print("This surface combination is not implemented:", ltype, surf_subtype, squarelens)
            return None, None, None, None, None
    else:
        if ltype == 'flat' or surf_subtype == 'flat': #flat surface case
            v, f, n = add_flat_surface(lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape)
            v_outline = v # trivial
        elif ltype == 'rotational' and surf_subtype == 'spherical':
            v, f, n = add_spherical_surface(srad, lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape, lrad_ext=lrad_ext)
            v_outline = v[-N2:]
        elif ltype == 'rotational' and surf_subtype != 'spherical':
            v, f, n = add_sagsurface_circular(srad, k, A, lrad, N1, N2, surf_rotation=surf_rotation, zadd=zadd, nVerts=nVerts, surftype='aspheric', dshape=dshape, lrad_ext=lrad_ext)
            v_outline = v[-N2:]
        elif ltype == 'cylindrical' and surf_subtype == 'spherical':
            v, f, n = add_sagsurface_circular(srad, k, A, lrad, N1, N2, surf_rotation=surf_rotation, zadd=zadd, nVerts=nVerts, surftype='cylindrical', dshape=dshape, lrad_ext=lrad_ext)
            v_outline = v[-N2:]
        elif ltype == 'cylindrical' and surf_subtype != 'spherical':
            v, f, n = add_sagsurface_circular(srad, k, A, lrad, N1, N2, surf_rotation=surf_rotation, zadd=zadd, nVerts=nVerts, surftype='acylindrical', dshape=dshape, lrad_ext=lrad_ext)
            v_outline = v[-N2:]
            """
        elif ltype == 'cylindricX' and surf_subtype == 'spherical':
            v, f, n = sfc.add_spherical_surface(srad, lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape, lrad_ext=lrad_ext, cylinderaxis='X')
        elif ltype == 'cylindricY' and surf_subtype == 'spherical':
            v, f, n = sfc.add_spherical_surface(srad, lrad, N1, N2, zadd=zadd, nVerts=nVerts, dshape=dshape, lrad_ext=lrad_ext, cylinderaxis='Y')
            """
        elif ltype == 'toric' and surf_subtype == 'spherical':
            v, f, n = add_sagsurface_circular(srad, k, A, lrad, N1, N2, zadd=zadd, nVerts=nVerts,
                    rad2=RY2, k2=None, A2=None,
                    surf_rotation=surf_rotation, surftype='toric',dshape=dshape, lrad_ext=lrad_ext)
            v_outline = v[-N2:]
        else:
            print("This surface combination is not implemented:", ltype, surf_subtype, squarelens)
            return None, None, None, None, None

    N_inside_sq = 'this_will_fail' # variable should later be overwritten, this way i will catch errors more easily
    if squarelens:
        # this variable is used later for side face generation
        if ltype == 'flat' or surf_subtype == 'flat':
            N_inside_sq = 2
        else:
            N_inside_sq = N2

    # reorder the coordinate system to Blender covnention
    v = [[t[2], t[0], t[1]] for t in v]
    n = [[t[2], t[0], t[1]] for t in n]
    v_outline = [[t[2], t[0], t[1]] for t in v_outline]

    return v, f, n, N_inside_sq, v_outline