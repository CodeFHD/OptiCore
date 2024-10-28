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

import numpy as np

def check_surface(r, lrad, flrad, k, A):
    ssig = 1 - 2*(r < 0) # this ensures ssig == 1 for r == 0
    
    # determine the type of surface w.r.t. the intersection algorithms:
    hasrad = r != 0
    hasconic = k != 0 and k is not None
    haspoly = not (np.all(np.array(A) == 0) or np.all(np.array(A) == None) or A is None)
    if not hasrad and not haspoly:
        surfshape = 'flat'
    elif not hasrad and haspoly:
        surfshape = 'polynomic'
    elif not hasconic and not haspoly:
        surfshape = 'spherical'
    elif hasrad and not haspoly and k == -1:
        surfshape = 'parabolic'
    elif hasconic and not haspoly:
        surfshape = 'conic'
    else: 
        surfshape = 'aspheric'

    # get parameters depending on surface type
    if surfshape == 'flat':
        # no flange needed for a flat surface
        hasfl = 0
        lrad_surf = lrad
    else:
        hasfl = flrad > 0.01*lrad
        if flrad > 0.99*lrad: flrad = 0.99*lrad
        lrad_surf = lrad - flrad if hasfl else lrad
        if k <= -1:
            # in this case 1+k is smaller than 1 and the situation is safe anyways
            pass
        elif lrad_surf > abs(r)/np.sqrt(1+k):
            lrad_surf = 0.9999*abs(r)/np.sqrt(1+k)
            flrad = lrad - lrad_surf
            hasfl = 1
            
    return lrad_surf, hasfl, flrad, ssig