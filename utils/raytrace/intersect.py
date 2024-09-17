import numpy as np

# Ray-Sphere intersection
def sphere_intersect(O, D, C, r, direction=1):
    """
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
    #calc intersect paramters
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
    #print('N: ', N[16])
    #print('P: ', P[16])
    return P, N
    
# Combine sphere_intersect with a check for aperture to create lens intersection
def lens_intersect(O, D, C, r, rad, direction=1):#, surface=-1):
    P, N = sphere_intersect(O, D, C, r, direction=direction)#, surface=surface)
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
def triangle_intersect(O, D, Tri):
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
    P = O + (t*D.T).T

    # Check if points are within triangle
    e0, e1, e2 = Tri[1] - Tri[0], Tri[2] - Tri[1], Tri[0] - Tri[2]
    p0, p1, p2 = P - Tri[0], P - Tri[1], P - Tri[2]
    s0 = np.sign(np.einsum('ij,ij->i', N, np.cross(e0, p0))).astype(int)
    s1 = np.sign(np.einsum('ij,ij->i', N, np.cross(e1, p1))).astype(int)
    s2 = np.sign(np.einsum('ij,ij->i', N, np.cross(e2, p2))).astype(int)
    #print(s0, s1, s2)
    #s = np.logical_and(s0 == s1, s0 == s2)
    s = s0*s1*s2
    s = s > -0.5

    P[~s] = float('nan')

    idx_fail = np.isnan(P[:,0])

    return P, N, idx_fail

def rectangle_intersect(O, D, Quad):
    """
    This function splits the rectangle into two Tris,
    checks intersection with each,
    and combines the result.
    Quad is shape (4, 3), and should be ordered
    """
    Tri1 = np.array([Quad[0], Quad[1], Quad[3]])
    Tri2 = np.array([Quad[1], Quad[2], Quad[3]])
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
    according to Snell's law for refraction.
    """
    NdotD = np.einsum('ij,ij->i', N,D)
    s = np.sign(-NdotD).astype(int)

    # get an angle difference accroding to Snell's law
    theta1 = np.arccos(np.abs(NdotD))
    theta2 = np.arcsin(np.sin(theta1)*n1/n2)
    dtheta = (theta2 - theta1)*s
    # For testing: print incoming and target outgoing angle
    # print(f'{n1:.4f} {n2:.4f} {theta1[0]*180/np.pi:.4f} {theta2[0]*180/np.pi:.4f}')
    
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
    # For testing: print outgoing angle as computed. Should be equal to target from above. (minus floating point errors)
    # print(f'{np.arccos(NdotD[0])*180/np.pi:.4f}')
        
    return Dnew

# Reflection
def reflect_ray(D, N):
    NdotD = np.einsum('ij,ij->i', N, D)
    s = np.sign(-NdotD).astype(int)
    Dnew = (D.T - N.T*2*NdotD).T
    return Dnew