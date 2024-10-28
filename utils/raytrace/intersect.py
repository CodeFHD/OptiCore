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

def intersect_sphere(O, D, C, r, k):
    n_rays = O.shape[0]
    ones = np.ones(n_rays)
    
    # Step 1: transform ray into local coordinate system 
    O = O - C
    
    # Step 2: intersect
    Ox, Oy, Oz = O[:,0], O[:,1], O[:,2]
    Dx, Dy, Dz = D[:,0], D[:,1], D[:,2]
    a = np.einsum('ij,ij->i', D, D)/r
    b = 2.*(np.einsum('ij,ij->i', O, D)/r - Dz)
    c = np.einsum('ij,ij->i', O, O)/r - 2*Oz
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
    Ox, Oy, Oz = O[:,0], O[:,1], O[:,2]
    Dx, Dy, Dz = D[:,0], D[:,1], D[:,2]
    a = 1./r*(Dx*Dx + Dy*Dy + (1+k)*Dz*Dz)
    b = 2.*(1./r*(Ox*Dx + Oy*Dy + (1+k)*Oz*Dz) - Dz)
    c = 1./r*(Ox*Ox + Oy*Oy + (1+k)*Oz*Oz) - 2*Oz
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
    dz = 1./r*(2*N1*N2 + N0)/(N1*N2*N2)

    N = np.column_stack((-P[:,0]*dz, -P[:,1]*dz, ones))
    N = (N.T/np.sqrt(np.einsum('ij,ij->i',N,N))).T # normalize the vector
    
    # Step 4: Transform point and normal back into original coordiante system
    P = P + C
        
    return P, N
    
# Combine sphere_intersect with a check for aperture to create lens intersection
def lens_intersect(O, D, C, r, rad, k=0, A=[0], surfshape='spherical', direction=1):#, surface=-1):
    
    if surfshape == 'flat':
        P , N, idx_fail = circle_intersect(O, D, C[2], rad)
    elif surfshape == 'conic':
        P, N = intersect_conic(O, D, C, r, k)
    else:
        P, N = intersect_sphere(O, D, C, r, 0)
        #P, N = sphere_intersect(O, D, C, r, direction=direction)#, surface=surface)
    #check if within lens radius, else set NaN
    diff1 = P[:,0]-C[0]
    diff2 = P[:,1]-C[1]
    Prad = np.sqrt(diff1*diff1 + diff2*diff2)
    #TODO: this could alternatively be swapped for marking only in idx_fail
    P[Prad > rad] = float('nan')
    idx_fail = np.isnan(P[:,0])
    return P, N, idx_fail

# Circle intersect is used e.g. for planar faces
def circle_intersect(O, D, zcirc, r, pass_inside=True):
    """
    Intersection with circular face
    oriented perpendicular to z-axis
    
    O is shape (N,3)
    D is shape (N,3)
    zcirc is float
    r is float  
    """
    #get Ray-t to circle assuming it is a vertical plane at z-location zdet
    t = (zcirc - O[:,2])/D[:,2]
    
    #calc intersection points
    P = O + (t*D.T).T
    
    #check if (not) within boundaries
    if pass_inside:
        idx_fail = np.einsum('ij,ij->i',P[:,:2],P[:,:2]) > r**2
        P[idx_fail] = float('nan')
    else:
        idx_fail = np.einsum('ij,ij->i',P[:,:2],P[:,:2]) < r**2
        P[idx_fail] = float('nan')

    #surface normal is in negative direction by definition
    N = np.array([[0,0,-1.] for i in range(P.shape[0])])
    
    return P, N, idx_fail

# Ray-Triangle intersection
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
    e0, e1, e2 = Tri[1] - Tri[0], Tri[2] - Tri[1], Tri[0] - Tri[2]
    p0, p1, p2 = P - Tri[0], P - Tri[1], P - Tri[2]
    s0 = np.sign(np.einsum('ij,ij->i', N, np.cross(e0, p0))).astype(int)
    s1 = np.sign(np.einsum('ij,ij->i', N, np.cross(e1, p1))).astype(int)
    s2 = np.sign(np.einsum('ij,ij->i', N, np.cross(e2, p2))).astype(int)
    s = s0*s1*s2
    s = s > -0.5

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

def aperture_intersect(O, D, r_ap, zap, n_blades=7):
    if n_blades < 3:
        #print('WARNING: Number of Blades below 3 doesnt make sense! Tracing circular aperture')
        n_blades = 0
    if n_blades == 0:
        P, N, idx_fail = circle_intersect(O, D, zap, r_ap, pass_inside=False)
        return P, N, idx_fail
        
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
    
    idx_fail = dist <= r_ap
    P[idx_fail] = float('nan')

    #surface normal is in negative direction by definition
    N = np.array([[0,0,-1.] for i in range(P.shape[0])])
    
    return P, N, idx_fail


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