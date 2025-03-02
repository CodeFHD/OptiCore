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

# Ray-Sphere intersection
def sphere_intersect(O, D, C, r, direction=1):
    """
    OBSOLETE AT THE MOMENT, KEEP FOR REFERENCE

    A good algorithm description can be found at:
    https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    
    O is shape (N,3) 
    D is shape (N,3)
    C is shape (3,)
    r is float
    direction is either +1 or -1
    surface is either +1 or -1
    """
    surface = np.sign(r).astype(int)
    # calc intersect paramters
    L = C - O
    absL2 = np.einsum('ij,ij->i', L, L) # Dot product
    tca = np.einsum('ij,ij->i', L, D) # Dot product
    d2 = absL2 - tca**2 # Pythagoras
    thc = np.sqrt(r**2 - d2) # Pythagoras
    if direction*surface == -1:
        t = np.abs(tca + thc)
    else:
        t = np.abs(tca - thc)
    td = (t*D.T).T # The transpose here is because of the shapes of the arrays
    P = O + td # intersection point
    N = P - C # surface normal at intersection point
    # Normalize the surface-normals
    N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T
    return P, N

def intersect_cylinder(O, D, C, surf_rotation, r):
    """
    C is a point along the cylinder axis
    surf_rotation is the angle the cylinder is rotated around the optical axis
    b is the unit vector along the cylinder axis (built from surf_rotation)
    r is the cylunder radius
    # tb is the length along c so a point M which marks the circle center with the intersection point P.
    
    https://stackoverflow.com/questions/73866852/ray-cylinder-intersection-formula
    """
    n_rays = O.shape[0]
    ones = np.ones(n_rays)
    # Explanation of surfrot: at 0 deg rotation the curvature is along the x-axis (y-axis displayed in Blender)
    # The axis perpendicular to this, i.e. along the y-axis
    # surfrot is defined to rotate counterclockwise , ending 
    b_cyl = [np.sin(-surf_rotation), np.cos(-surf_rotation), 0]
    b_cyl = np.array(b_cyl)
    
    # Step 1: transform ray into local coordinate system
    O = O - C

    # Step 2: intersect
    Ox, Oy, Oz = O[:,0], O[:,1], O[:,2]
    Dx, Dy, Dz = D[:,0], D[:,1], D[:,2]
    a = np.einsum('ij,ij->i', D, D) - np.einsum('ij,j->i', D, b_cyl)**2
    b = 2*(np.einsum('ij,ij->i', D, O) - np.einsum('ij,j->i', D, b_cyl)*np.einsum('ij,j->i', O, b_cyl) -Dz*r)
    c = np.einsum('ij,ij->i', O, O) - np.einsum('ij,j->i', O, b_cyl)**2- 2*Oz*r
    root2 = b*b - 4*a*c
    # get the two possible roots
    t1 = (-b - np.sqrt(root2)) / (2*a)
    t2 = (-b + np.sqrt(root2)) / (2*a)
    # special case of a == 0 breaks quadratic roots
    if np.any(a==0):
        idxa = a==0
        ta = -c/b
        t1[idxa] = ta[idxa]

    # test both t-values and determine plausible intersection
    # intersection closer to Origin (0, 0, 0) will win here because we are working in local coordinates    
    # Get Points
    td1 = (t1*D.T).T # The transpose here is because of the shapes of the arrays
    P = O + td1 # intersection point
    td2 = (t2*D.T).T # The transpose here is because of the shapes of the arrays
    P2 = O + td2 # intersection point
    # Distance to origin
    absP1 = np.einsum('ij,ij->i',P,P)
    absP2 = np.einsum('ij,ij->i',P2,P2)
    idx = absP2 < absP1
    P[idx] = P2[idx]
    t1[idx] = t2[idx]

    # Step 3: Calculate Normal
    m = np.einsum('ij,j->i', D, b_cyl)*t1 + np.einsum('ij,j->i', O, b_cyl)
    C2 = [0, 0, r] + (m*b_cyl[:,np.newaxis]).T
    # C2 = [0, 0, r]
    N = P - C2 # surface normal at intersection point
    # Normalize the surface-normals
    N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T

    # Step 4: Transform point and normal back into original coordiante system
    P = P + C
        
    return P, N

def intersect_sphere(O, D, C, r, k):
    n_rays = O.shape[0]
    ones = np.ones(n_rays)
    
    # Step 1: transform ray into local coordinate system 
    O = O - C
    
    # Step 2: intersect
    Oz = O[:,2]
    Dz = D[:,2]
    a = 1 # np.einsum('ij,ij->i', D, D) # This wokrs only because/when D is normalized
    b = 2.*(np.einsum('ij,ij->i', O, D) - Dz*r)
    c = np.einsum('ij,ij->i', O, O) - 2*Oz*r
    root2 = b*b - 4*c # 4*a*c # remember a == 1
    # get the two possible roots
    t1 = (-b - np.sqrt(root2)) / 2 # (2*a)
    t2 = (-b + np.sqrt(root2)) / 2 # (2*a)
    # special case of a == 0 breaks quadratic roots
    # doesn't apply here because a can't be zero
    """
    if np.any(a==0):
        idxa = a==0
        ta = -c/b
        t1[idxa] = ta[idxa]
    """
        
    # test both t-values and determine plausible intersection
    # intersection closer to Origin (0, 0, 0) will win here because we are working in local coordinates    
    # Get Points
    td1 = (t1*D.T).T # The transpose here is because of the shapes of the arrays
    P = O + td1 # intersection point
    td2 = (t2*D.T).T # The transpose here is because of the shapes of the arrays
    P2 = O + td2 # intersection point
    # Distance to origin
    absP1 = np.einsum('ij,ij->i',P,P)
    absP2 = np.einsum('ij,ij->i',P2,P2)
    idx = absP2 < absP1
    P[idx] = P2[idx]
    
    # Step 3: Calculate Normal
    C2 = [0, 0, r]
    N = P - C2 # surface normal at intersection point
    # Normalize the surface-normals
    N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T
    
    # Step 4: Transform point and normal back into original coordiante system
    P = P + C
        
    return P, N

def intersect_conic(O, D, C, r, k):
    n_rays = O.shape[0]
    ones = np.ones(n_rays)
    
    # Step 1: transform ray into local coordinate system 
    O = O - C
    
    # Step 2: intersect
    Oz = O[:,2]
    Dz = D[:,2]
    a = 1 + k*Dz*Dz
    b = 2.*(np.einsum('ij,ij->i', O, D) - Dz*r + k*Oz*Dz)
    c = np.einsum('ij,ij->i', O, O) - 2*Oz*r + k*Oz*Oz
    root2 = b*b - 4*a*c
    # get the two possible roots
    t1 = (-b - np.sqrt(root2)) / (2*a)
    t2 = (-b + np.sqrt(root2)) / (2*a)
    # special case of a == 0 breaks quadratic roots
    if np.any(a==0):
        idxa = a==0
        ta = -c/b
        t1[idxa] = ta[idxa]
    
    # test both t-values and determine plausible intersection
    # intersection closer to Origin (0, 0, 0) will win here because we are working in local coordinates    
    # Get Points
    td1 = (t1*D.T).T # The transpose here is because of the shapes of the arrays
    P = O + td1 # intersection point
    td2 = (t2*D.T).T # The transpose here is because of the shapes of the arrays
    P2 = O + td2 # intersection point
    # Distance to origin
    absP1 = np.einsum('ij,ij->i',P,P)
    absP2 = np.einsum('ij,ij->i',P2,P2)
    idx = absP2 < absP1
    P[idx] = P2[idx]

    # Step 3: Calculate Normal
    R2 = P[:,0]*P[:,0] + P[:,1]*P[:,1]
    N0 = (1+k)*R2/r/r
    N1 = np.sqrt(1-N0)
    N2 = 1+N1
    dzdr = 1./r*(2*N1*N2 + N0)/(N1*N2*N2)

    N = np.column_stack((-P[:,0]*dzdr, -P[:,1]*dzdr, ones))
    N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T # normalize the vector
    
    # Step 4: Transform point and normal back into original coordiante system
    P = P + C
        
    return P, N
    
# Combine sphere_intersect with a check for aperture to create lens intersection
def lens_intersect(O, D, C, r, rad, k=0, A=[0], surf_rotation=0, surfshape='spherical', direction=1, ):#, surface=-1):
    if surfshape == 'flat':
        P, _, N, idx_fail, _ = circle_intersect(O, D, C[2], rad)
    elif surfshape == 'conical':
        P, N = intersect_conic(O, D, C, r, k)
    elif surfshape == 'cylindrical':
        P, N = intersect_cylinder(O, D, C, surf_rotation, r)
    else:
        P, N = intersect_sphere(O, D, C, r, 0)
        #P, N = sphere_intersect(O, D, C, r, direction=direction)#, surface=surface)
    #check if within lens radius, else set NaN
    diff1 = P[:,0]-C[0]
    diff2 = P[:,1]-C[1]
    Prad = np.sqrt(diff1*diff1 + diff2*diff2)
    P[Prad > rad] = float('nan') # TODO: remove this, only create idx_fail, and let rays.update() take care of the rest ?
    idx_fail = np.isnan(P[:,0])
    return P, N, idx_fail

# Circle intersect is used e.g. for planar faces
def circle_intersect(O, D, zcirc, r, pass_inside=True, pass_outside=False):
    """
    Intersection with circular face
    oriented perpendicular to z-axis
    
    O is shape (N,3)
    D is shape (N,3)
    zcirc is float
    r is float  
    """
    if not pass_inside and not pass_outside:
        print("[OC] WARNING: circle_intersect() not specified with pass_inside or pass_outside. Defaulting to pass_inside.")
        pass_inside = True
    P_inside = None
    P_outside = None

    # get Ray-t to circle assuming it is a vertical plane at z-location zdet
    t = (zcirc - O[:,2])/D[:,2]
    
    # calc intersection points
    P = O + (t*D.T).T
    
    # check if (not) within boundaries
    if pass_inside:
        P_inside = np.array(P)
        idx_fail_i = np.einsum('ij,ij->i',P[:,:2],P[:,:2]) > r**2
        P_inside[idx_fail_i] = float('nan')
    if pass_outside:
        P_outside = np.array(P)
        idx_fail_o = np.einsum('ij,ij->i',P[:,:2],P[:,:2]) <= r**2
        P_outside[idx_fail_o] = float('nan')

    # surface normal is in negative direction by definition
    N = np.array([[0, 0, -1.] for i in range(P.shape[0])])

    if P_inside is None: idx_fail_i = None
    if P_outside is None: idx_fail_o = None
    return P_inside, P_outside, N, idx_fail_i, idx_fail_o

# Ray-Triangle intersection
def is_point_in_tri(P, Tri, N):
    e0, e1, e2 = Tri[1] - Tri[0], Tri[2] - Tri[1], Tri[0] - Tri[2]
    p0, p1, p2 = P - Tri[0], P - Tri[1], P - Tri[2]
    s0 = np.sign(np.einsum('j,ij->i', N, np.cross(e0, p0))).astype(int) > -0.5
    s1 = np.sign(np.einsum('j,ij->i', N, np.cross(e1, p1))).astype(int) > -0.5
    s2 = np.sign(np.einsum('j,ij->i', N, np.cross(e2, p2))).astype(int) > -0.5
    s = s0 & s1 & s2
    return s

def triangle_intersect(O, D, Tri, return_t=False):
    """
    O is shape (N,3) 
    D is shape (N,3)
    Tri is shape (3,3), each row being one vertex
    Triangle vertices should be ordered so that the surface normal is alogn the negative z-direction,
    i.e. right-handed
    """
    # compute the triangles surface normal
    n1 = Tri[1] - Tri[0]
    n2 = Tri[2] - Tri[0]
    n = np.cross(n1, n2)
    n = n/np.sqrt(np.dot(n,n))
    # create array the same size of input rays
    N = np.array([n for i in range(O.shape[0])])

    # d comes from the 4-coefficeint equation of a plane
    d = -1*np.dot(n, Tri[0])

    # compute ray-plane intersection
    t = -1 * (np.einsum('ij,j->i', O, n) + d) / np.einsum('ij,j->i', D, n)
    # filter huge values that occur when the ray is parallel to the plane
    # but somehow the check within triangle below returns True
    t[t > 1e9] = float('nan')
    P = O + (t*D.T).T

    # Check if points are within triangle
    s = is_point_in_tri(P, Tri, n)

    P[~s] = float('nan')

    idx_fail = np.isnan(P[:,0])
    
    if return_t:
        t[idx_fail] = float('nan')
        return t
    return P, N, idx_fail

def rectangle_intersect(O, D, Quad, return_t=False):
    """
    This function splits the rectangle into two Tris,
    checks intersection with each,
    and combines the result.
    Quad is shape (4, 3), and should be properly ordered
    """
    if np.all(np.isnan(O)):
        return O, D, np.isnan(O[:,0])
    Tri1 = np.array([Quad[0], Quad[1], Quad[3]])
    Tri2 = np.array([Quad[1], Quad[2], Quad[3]])
    if return_t:
        t1 = triangle_intersect(O, D, Tri1, return_t=True)
        t2 = triangle_intersect(O, D, Tri2, return_t=True)
        idx = np.where(np.isnan(t1))[0]
        t1[idx] = t2[idx]
        return t1
    P1, N1, idx_fail1 = triangle_intersect(O, D, Tri1)
    P2, N2, idx_fail2 = triangle_intersect(O, D, Tri2)
    P1[~idx_fail2] = P2[~idx_fail2]
    N1[~idx_fail2] = N2[~idx_fail2]
    idx_fail1[~idx_fail2] = idx_fail2[~idx_fail2]
    return P1, N1, idx_fail1

def aperture_intersect(O, D, r_ap, zap, n_blades=7, pass_inside=False):
    if n_blades < 3:
        #print('WARNING: Number of Blades below 3 doesnt make sense! Tracing circular aperture')
        n_blades = 0
    if n_blades == 0:
        P_inside, P_outside, N, idx_fail_i, idx_fail_o = circle_intersect(O, D, zap, r_ap, pass_inside=pass_inside, pass_outside=True)
        return P_inside, P_outside, N, idx_fail_i, idx_fail_o
    
    P_inside = None
    P_outside = None
        
    #get Ray-t to plane assuming it is a vertical plane at z-location zap
    t = (zap - O[:,2])/D[:,2]
    
    #calc intersection points
    P = O + (t*D.T).T
    
    amax = 2*np.pi/n_blades
    a = np.arctan2(P[:,1], P[:,0])
    section = (a + amax/2)//amax
    section = section*amax
    
    bladedir = np.vstack((np.sin(section), np.cos(section))).T
    
    dist = np.einsum('ij,ij->i', P[:,:2], bladedir)
    
    if pass_inside:
        P_inside = np.array(P)
        idx_fail_i = dist > r_ap
        P_inside[idx_fail_i] = float('nan')
    # pass_outside always True in this function
    P_outside = np.array(P)
    idx_fail_o = dist <= r_ap
    P_outside[idx_fail_o] = float('nan')

    #surface normal is in negative direction by definition
    N = np.array([[0,0,-1.] for i in range(P.shape[0])])
    
    return P_inside, P_outside, N, idx_fail_i, idx_fail_o


#################################################

# Refraction
def refract_ray(D, N, n1, n2):
    """
    This function uses Rodrigues' rotation formula
    in accordance with Snell's law for refraction.
    """
    
    if np.all(np.isnan(D)):
        return D
    
    # print('refraction with', n1, n2)

    NdotD = np.einsum('ij,ij->i', N,D)
    s = np.sign(-NdotD).astype(int)

    # get an angle difference accroding to Snell's law
    theta1 = np.arccos(np.abs(NdotD))
    theta2 = np.arcsin(np.sin(theta1)*n1/n2)
    dtheta = (theta2 - theta1)*s
    
    k = -1*np.cross(N, D)
    absk = np.sqrt(np.einsum('ij,ij->i',k,k))
    ak0 = absk == 0 #the computation above fails if ray is along the normal.
    k[ak0] = [0,0,1]
    absk[ak0] = 1
    k = k.T/absk
    k = k.T

    kcrossD = np.cross(k,D)
    kdotD = np.einsum('ij,ij->i', k,D)
    cdtheta = np.cos(dtheta)
    sdtheta = np.sin(dtheta)
    kcvsdt = kcrossD.T*sdtheta
    Dcdt = D.T*cdtheta
    kkdotD = k.T*kdotD
    kkdotDcdt = kkdotD*(1-cdtheta)
    
    Dnew = 1*Dcdt.T + kcvsdt.T + kkdotDcdt.T
    NdotD = np.abs(np.einsum('ij,ij->i', -N,Dnew))
        
    return Dnew

# Reflection
def reflect_ray(D, N):
    if np.all(np.isnan(D)):
        return D
    NdotD = np.einsum('ij,ij->i', N, D)
    s = np.sign(-NdotD).astype(int)
    Dnew = (D.T - N.T*2*NdotD).T
    return Dnew