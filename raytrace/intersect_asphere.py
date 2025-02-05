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

from .intersect import rectangle_intersect
from ..surface.radialprofiles import get_z_evenasphere, get_N_evenasphere
from ..surface.toric import get_z_toric, get_N_toric

"""
PARKED CODE for spossible use of Newton-Raphson in the future

def dzdt(k, A, O, D, t):
    pass

def dzdr(c, k, A, r):
    pass

def dz1dr(c, k, r):
    pass

def dz2dr(A, r):
    r_dz2dr = sum([2*(i+2)*A[i]*r**(2*(i+2) - 1) for i in range(len(A))])
    return r_dz2dr

def dr2dt(O, D, t):
    # t is shape (N)
    # O,D are shape (N,3)
    
    D2 = D*D
    OD = O*D
    r_dr2dt = 2*(t*(D2[:,0] + D2[:,1]) + (OD[:,0] + OD[:,1]))
    return r_dr2dt

def drdt(O, D, t, r):
    D2 = D*D
    OD = O*D
    O2 = O*O
    D2s = (D2[:,0] + D2[:,1])
    ODs = (OD[:,0] + OD[:,1])
    O2s = (O2[:,0] + O2[:,1])
    r2 = t*t*D2s + 2*t*ODs + O2s
    r = np.sqrt(r2)
    r_drdt = (t*D2s + ODs) / r
    return r_drdt
"""    

# configuration for the algorithms
N1 = 17
N2 = 16 # Not needed as long as rotational symmetry of lens is assumed
N_BISECTION = 12
ACCURACY = 1e-5
MAX_ITERATIONS = 6

def _get_bbox(lrad, zmin, zmax):
    # Quads for the 6 box sides
    V1 = [-lrad, -lrad, zmin]
    V2 = [lrad, -lrad, zmin]
    V3 = [lrad, lrad, zmin]
    V4 = [-lrad, lrad, zmin]
    V5 = [-lrad, -lrad, zmax]
    V6 = [lrad, -lrad, zmax]
    V7 = [lrad, lrad, zmax]
    V8 = [-lrad, lrad, zmax]
    # generate the sides Faces
    Qfront = [V1, V2, V3, V4]
    Qback = [V5, V6, V7, V8]
    Qleft = [V1, V4, V8, V5]
    Qright = [V2, V3, V7, V6]
    Qtop = [V4, V3, V7, V8]
    Qbottom = [V1, V2, V6, V5]

    return Qfront, Qback, Qleft, Qright, Qtop, Qbottom

def intersect_asphere(O, D, C, lrad, rad, k, A):
    n_rays = O.shape[0]
    ones = np.ones(n_rays)
    
    # Step 1: transform ray into local coordinate system 
    O = O - C
    
    ###########################################################################
    
    # generate bounding box
    r = np.linspace(-lrad, lrad, N1, endpoint = True)
    z = get_z_evenasphere(r**2, rad, k, A)
    zmin = 1.2*min(z) # 20 percent margin because of the coarse sampling of the surface
    zmax = 1.2*max(z)
    Qfront, Qback, Qleft, Qright, Qtop, Qbottom = _get_bbox(lrad, zmin, zmax)
    
    # intersect ray with bbox to find bracket for root-finding
    # initialize min and max t values, nan by defaul
    t0 = np.full((n_rays), float('nan'))
    t1 = np.full((n_rays), float('nan'))
    # try front-back first. should account for most rays
    t_front = rectangle_intersect(O, D, Qfront, return_t=True)
    t_back = rectangle_intersect(O, D, Qback, return_t=True)
    idx_front = np.isnan(t_front)
    idx_back = np.isnan(t_back)
    # ASSUMPTION: front can only be first, back only last == no rays travelling backwards
    # TODO: Deal with Ghosts, mirrors, and catadioptrics later after you confirmed this works
    t0[~idx_front] = t_front[~idx_front]
    t1[~idx_back] = t_back[~idx_back]
    # check sides for any where front and/or back failed
    if np.any(np.isnan(t0)) or np.any(np.isnan(t1)):
        t_left = rectangle_intersect(O, D, Qleft, return_t=True)
        t_right = rectangle_intersect(O, D, Qright, return_t=True)
        t_top = rectangle_intersect(O, D, Qtop, return_t=True)
        t_bottom = rectangle_intersect(O, D, Qbottom, return_t=True)
        idx_left = np.isnan(t_left)
        idx_right = np.isnan(t_right)
        idx_top = np.isnan(t_top)
        idx_bottom = np.isnan(t_bottom)
        # compare t to t_front and t_back to determin if hit must be a first or second
        t_all = np.column_stack((t_left, t_right, t_top, t_bottom))
        t_min = np.nanmin(t_all, axis=1)
        t_max = np.nanmax(t_all, axis=1)
        t0[idx_front] = t_min[idx_front]
        t1[idx_back] = t_max[idx_back]
    
    # use iterative bisection method to find t    
    n_iteration = 0
    while n_iteration < MAX_ITERATIONS:
        n_iteration = n_iteration + 1
        # generate the test points
        t_guess = np.linspace(t0, t1, N_BISECTION, endpoint=True)
        DT = np.einsum('ij,ki->ijk', D, t_guess)
        P_guess = O[:,:,np.newaxis] + DT
        # evaluate the z-distance from the surface
        r2_guess = np.einsum('ijk,ijk->ik', P_guess[:,:2,:], P_guess[:,:2,:])
        z_guess = P_guess[:,2,:]
        z_surf = get_z_evenasphere(r2_guess, rad, k, A)
        dz_guess = z_guess - z_surf
        dz_guess_abs = np.abs(dz_guess)
        dz_min = np.min(dz_guess_abs, axis=1)#[:,np.newaxis]
        if not np.any(dz_min > ACCURACY):
            break
        # update guess for t
        # idx_min = np.argmin(dz_guess_abs, axis = 1)
        # t = np.choose(idx_min, t_guess)
        idx_gt = (dz_guess > 0).argmax(1)
        idx_gt[idx_gt == 0] = N_BISECTION - 1
        t0 = np.choose(idx_gt-1, t_guess)# - ACCURACY
        t1 = np.choose(idx_gt, t_guess)# + ACCURACY
    
    # final evaluation
    idx_min = np.argmin(dz_guess, axis = 1)
    t = np.choose(idx_min, t_guess)
    # discard rays that did not meet ACCURACY
    idx_accuracy = dz_min > ACCURACY
    t[idx_accuracy] = float('nan')

    # nfail = sum(idx_accuracy)
    # print(f'Aspheric trace finished after {n_iteration}/{MAX_ITERATIONS} iterations with {nfail}/{n_rays} failures')
    
    td = (t*D.T).T # The transpose here is because of the shapes of the arrays
    P = O + td # intersection point
    # C2 = [0, 0, rad]
    
    # calculate surface normal
    r = np.sqrt(np.einsum('ij,ij->i', P[:,:2], P[:,:2]))
    b = np.arctan2(P[:,1], P[:,0])
    N = get_N_evenasphere(r, b, rad, k, A, returnformat='seqrt')
    N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T
    # Flip N towards incoming ray direction
    s = np.sign(np.einsum('ij,ij->i', N, D))
    s[s==0] = 1
    N = (N.T*s).T
    
    # TODO optional: use bisection method to get initial guess, then use Newton-Raphson for refinement

    ###########################################################################
    
    # Step 4: Transform point and normal back into original coordiante system
    P = P + C
    
    # Step 5: Check if Point lies within defined aperture
    diff1 = P[:,0]-C[0]
    diff2 = P[:,1]-C[1]
    Prad = np.sqrt(diff1*diff1 + diff2*diff2)
    #TODO: this could alternatively be swapped for marking only in idx_fail
    P[Prad > lrad] = float('nan')

    
    idx_fail = np.isnan(P[:,0])
    return P, N, idx_fail 



def intersect_implicit(O, D, C, lrad, z_fun, z_fun_params, N_fun=None):
    """
    Intersection for a general implicit surface,
    i.e. one that takes the form f(x,y,z) = 0.
    This function needs to be supplied with a function reference (plus ccoefficients)
    z_fun so that z_fun(x, y, *z_fun_params) = z, 
    as well as a function N_fun that calculates the local surface normal.
    The parameter list z_fun_params shall also apply to N_fun.
    TODO: Add option to numerically obtain N_fun from z_fun when there is no (known) analytic solution.
    """
    n_rays = O.shape[0]
    ones = np.ones(n_rays)

    # Step 1: transform ray into local coordinate system 
    O = O - C

    ###########################################################################
    
    # generate bounding box
    x_test = [0]
    y_test = [0]
    for i in range(N1):
        r = lrad*i/(N1 - 1)
        phi_list = np.linspace(0, 2*np.pi, N2, endpoint=False)
        x_new = [r*np.cos(phi) for phi in phi_list]
        y_new = [r*np.sin(phi) for phi in phi_list]
        x_test = x_test + x_new
        y_test = y_test + y_new
    x_test, y_test = np.array(x_test), np.array(y_test)
    z_test = z_fun(x_test, y_test, *z_fun_params)
    zmin = 1.2*min(z_test) # 20 percent margin because of the coarse sampling of the surface
    zmax = 1.2*max(z_test)
    Qfront, Qback, Qleft, Qright, Qtop, Qbottom = _get_bbox(lrad, zmin, zmax)
    
    # intersect ray with bbox to find bracket for root-finding
    # initialize min and max t values, nan by defaul
    t0 = np.full((n_rays), float('nan'))
    t1 = np.full((n_rays), float('nan'))
    # try front-back first. should account for most rays
    t_front = rectangle_intersect(O, D, Qfront, return_t=True)
    t_back = rectangle_intersect(O, D, Qback, return_t=True)
    idx_front = np.isnan(t_front)
    idx_back = np.isnan(t_back)
    # ASSUMPTION: front can only be first, back only last == no rays travelling backwards
    # TODO: Deal with Ghosts, mirrors, and catadioptrics later after you confirmed this works
    t0[~idx_front] = t_front[~idx_front]
    t1[~idx_back] = t_back[~idx_back]
    # check sides for any where front and/or back failed
    if np.any(np.isnan(t0)) or np.any(np.isnan(t1)):
        t_left = rectangle_intersect(O, D, Qleft, return_t=True)
        t_right = rectangle_intersect(O, D, Qright, return_t=True)
        t_top = rectangle_intersect(O, D, Qtop, return_t=True)
        t_bottom = rectangle_intersect(O, D, Qbottom, return_t=True)
        idx_left = np.isnan(t_left)
        idx_right = np.isnan(t_right)
        idx_top = np.isnan(t_top)
        idx_bottom = np.isnan(t_bottom)
        # compare t to t_front and t_back to determin if hit must be a first or second
        t_all = np.column_stack((t_left, t_right, t_top, t_bottom))
        t_min = np.nanmin(t_all, axis=1)
        t_max = np.nanmax(t_all, axis=1)
        t0[idx_front] = t_min[idx_front]
        t1[idx_back] = t_max[idx_back]
    
    # use iterative bisection method to find t    
    n_iteration = 0
    while n_iteration < MAX_ITERATIONS:
        n_iteration = n_iteration + 1
        # generate the test points
        t_guess = np.linspace(t0, t1, N_BISECTION, endpoint=True)
        DT = np.einsum('ij,ki->ijk', D, t_guess)
        P_guess = O[:,:,np.newaxis] + DT
        # evaluate the z-distance from the surface
        #r2_guess = np.einsum('ijk,ijk->ik', P_guess[:,:2,:], P_guess[:,:2,:])
        x_guess = P_guess[:,0,:]
        y_guess = P_guess[:,1,:]
        z_guess = P_guess[:,2,:]
        z_surf = z_fun(x_guess, y_guess, *z_fun_params)
        #z_surf = get_z_evenasphere(r2_guess, rad, k, A)
        dz_guess = z_guess - z_surf
        dz_guess_abs = np.abs(dz_guess)
        dz_min = np.min(dz_guess_abs, axis=1)#[:,np.newaxis]
        if not np.any(dz_min > ACCURACY):
            break
        # update guess for t
        # idx_min = np.argmin(dz_guess_abs, axis = 1)
        # t = np.choose(idx_min, t_guess)
        idx_gt = (dz_guess > 0).argmax(1)
        idx_gt[idx_gt == 0] = N_BISECTION - 1
        t0 = np.choose(idx_gt-1, t_guess)# - ACCURACY
        t1 = np.choose(idx_gt, t_guess)# + ACCURACY
    
    # final evaluation
    idx_min = np.argmin(dz_guess, axis = 1)
    t = np.choose(idx_min, t_guess)
    # discard rays that did not meet ACCURACY
    idx_accuracy = dz_min > ACCURACY
    t[idx_accuracy] = float('nan')

    # nfail = sum(idx_accuracy)
    # print(f'Aspheric trace finished after {n_iteration}/{MAX_ITERATIONS} iterations with {nfail}/{n_rays} failures')
    
    td = (t*D.T).T # The transpose here is because of the shapes of the arrays
    P = O + td # intersection point
    
    # calculate surface normal
    x_final = P[:,0]
    y_final = P[:,1]
    N = N_fun(x_final, y_final, *z_fun_params, returnformat='seqrt')

    #r = np.sqrt(np.einsum('ij,ij->i', P[:,:2], P[:,:2]))
    #b = np.arctan2(P[:,1], P[:,0])
    #N = get_N_evenasphere(r, b, rad, k, A, returnformat='seqrt')
    # N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T
    # Flip N towards incoming ray direction
    s = np.sign(np.einsum('ij,ij->i', N, D))
    s[s==0] = 1
    N = (N.T*s).T
    
    # TODO optional: use bisection method to get initial guess, then use Newton-Raphson for refinement

    ###########################################################################
    
    # Step 4: Transform point and normal back into original coordiante system
    P = P + C
    
    # Step 5: Check if Point lies within defined aperture
    diff1 = P[:,0] # - C[0]
    diff2 = P[:,1] # - C[1]
    Prad = np.sqrt(diff1*diff1 + diff2*diff2)
    #TODO: this could alternatively be swapped for marking only in idx_fail
    P[Prad > lrad] = float('nan')
    
    idx_fail = np.isnan(P[:,0])
    return P, N, idx_fail 