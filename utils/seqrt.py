import bpy
import numpy as np

from mathutils import Vector

from bpy_extras.object_utils import AddObjectHelper, object_data_add

from . import rayfan

def rotate_rodrigues(D, N, n1, n2, refract=1, direction=1, surface=1):
    """
    Rodrigues rotation formula
    """
    NdotD = np.einsum('ij,ij->i', N,D)

    if surface == direction:
        theta1 = np.pi - np.arccos(NdotD)
    else:
        theta1 = np.arccos(NdotD)
       
    if refract == 1:
        theta2 = np.arcsin(np.sin(theta1)*n1/n2)
        dtheta = theta1 - theta2 
    else:
        dtheta = np.pi - 2*theta1
        dtheta *= surface

    #based on certain conditions, change the rotation angle
    #maybe this can be simplified    
    if refract == 1:
        if direction == 1:
            if surface == -1:
                dtheta *= -1
        else:
            if surface == 1:
                dtheta *= -1
    else:
        if direction == 1:
                dtheta*= -1
    
    #TODO: test for TIR
    k = np.cross(N, D)
    absk = np.sqrt(np.einsum('ij,ij->i',k,k))
    ak0 = absk == 0 #the computation fails if ray is along the normal. Treated at the end
    k[ak0] = [1,0,0]
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
    
    #correct direct incidence for reflection
    #Dnew *= refract
    if refract == -1:
        Dnew[ak0] = -1*Dnew[ak0]
        
    return Dnew

def plane_intersect(O, D, xdet):
    """
    Intersection with Detector
    
    O is shape (N,3)
    D is shape (N,3)
    xdet is float 
    """
    
    #get Ray-t to plane assuming it is a vertical plane at X-location xdet
    t = (xdet - O[:,0])/D[:,0]
    
    #calc intersection points
    P = O + (t*D.T).T
    
    return P

def sphere_intersect(O, D, C, r, rad, n1, n2, refract=1, direction=1):
    """
    A good algorithm description can be found at:
    https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    
    O is shape (N,3)
    D is shape (N,3)
    C is shape (3,)
    r, rad, n1, n2 are float
    refract and direction are either +1 or -1
    """
    
    surface = np.sign(r).astype(int)
    
    #t00 = time.perf_counter()
    
    #calc intersect paramters
    if True: #python/numpy implementation 
        L = C - O
        absL2 = np.einsum('ij,ij->i',L,L) #Dot product
        tca = np.einsum('ij,ij->i', L,D) #Dot product
        #TODO: here could be a check if tca < 0 return false (ray away from surface)
        d2 = absL2 - tca**2 #Pythagoras
        #TODO: here could be a check if d >= rad return false. Now automatically handeled by next line returning NaN
        thc = np.sqrt(r**2 - d2) #Pythagoras
        t = tca - direction*surface*thc
        #TODO: check if t <= 0 (?)
        td = D.T*t
        P = O + td.T
        N = P - C
    """
    else: #Using the C++ library
        P = np.zeros(O.shape,dtype = np.float32)
        N = np.zeros(O.shape,dtype = np.float32)
        suces = SeqRT.sphereintersect(np.ascontiguousarray(O, dtype=np.float32),
                                      np.ascontiguousarray(D, dtype=np.float32), 
                                      np.ascontiguousarray(C, dtype=np.float32), 
                                      r, rad, P, N)
    """
        
    #Normalize normals
    N = N.T/np.sqrt(np.einsum('ij,ij->i',N,N))
    N = N.T#*surface #*surface with C++ lib
    
    #if printdebug: print('Tracing part 1 took: ', time.perf_counter() - t00)
    
    #check if within lens radius, else set NaN
    diff1 = P[:,1]-C[1]
    diff2 = P[:,2]-C[2]
    Prad = np.sqrt(diff1*diff1 + diff2*diff2)
    P[Prad > rad] = float('nan')
    
    return P, N

def trace_rays(self, context):
    verts = []
    edges = []
    faces = []

    md = self.makedoublet
    nsurf = 2 + md

    i0 = 1
    i1 = self.ior
    if md:
        i2 = self.ior2
        i3 = 1
    else:
        i2 = 1

    lr = self.lensradius
    nrays = self.nrays
    if self.fantype == 'f2d':
        O, D = rayfan.rayfan2D(nrays, 0.9*lr, -1.0*lr)
    elif self.fantype == 'f3d':
        O, D = rayfan.rayfan3D(nrays, 0.9*lr, -1.0*lr)
        nrays = nrays*nrays - nrays + 1
    elif self.fantype == 'f3dr':
        O, D = rayfan.rayfan3D_uniformdiskrandom(nrays, 0.9*lr, -1.0*lr)
    elif self.fantype == 'f3dt':
        O, D = rayfan. rayfan3D_tri(nrays, 0.9*lr, -1.0*lr)
        nrays = nrays*nrays - nrays//2

    nVerts = 0

    #surface 1
    r = self.rad1
    surface = np.sign(r).astype(int)
    C = np.array([self.rad1, 0, 0])
    P, N = sphere_intersect(O, D, C, r, lr, i0, i1)
    D = rotate_rodrigues(D, N, i0, i1, 1, 1, surface)
    #make rays
    for i in range(nrays):
        verts.append(Vector(O[i]))
        verts.append(Vector(P[i]))
        edges.append([nVerts + 2*i, nVerts + 2*i + 1])
    #reinit
    nVerts = len(verts)
    O = np.array(P)

    #surface 2
    r = self.rad2
    surface = np.sign(r).astype(int)
    C = np.array([self.centerthickness + self.rad2, 0, 0])
    P, N = sphere_intersect(O, D, C, r, lr, i1, i2)
    D = rotate_rodrigues(D, N, i1, i2, 1, 1, surface)
    #make rays
    for i in range(nrays):
        verts.append(Vector(O[i]))
        verts.append(Vector(P[i]))
        edges.append([nVerts + 2*i, nVerts + 2*i + 1])
    #reinit
    nVerts = len(verts)
    O = np.array(P)

    if md:
        #surface 3
        r = self.rad3
        surface = np.sign(r).astype(int)
        C = np.array([self.centerthickness + self.centerthickness2 + self.rad3, 0, 0])
        P, N = sphere_intersect(O, D, C, r, lr, i2, i3)
        D = rotate_rodrigues(D, N, i2, i3, 1, 1, surface)
        #make rays
        for i in range(nrays):
            verts.append(Vector(O[i]))
            verts.append(Vector(P[i]))
            edges.append([nVerts + 2*i, nVerts + 2*i + 1])
        #reinit
        nVerts = len(verts)
        O = np.array(P)

    #detector
    P = plane_intersect(O, D, self.centerthickness + self.centerthickness2 + self.zdet)
    #make rays
    for i in range(nrays):
        verts.append(Vector(O[i]))
        verts.append(Vector(P[i]))
        edges.append([nVerts + 2*i, nVerts + 2*i + 1])
    #reinit
    nVerts = len(verts)
    O = np.array(P)

    #flip convention
    verts2 = []
    for v in verts:
        verts2.append(Vector((-1.*v[0],v[1],v[2])))

    #create mesh from verts and faces
    mesh = bpy.data.meshes.new(name="New Lens")
    mesh.from_pydata(verts2, edges, faces)
    obj = object_data_add(context, mesh, operator=self)