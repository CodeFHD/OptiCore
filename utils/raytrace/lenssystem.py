from . import glasscatalog

class Lenssystem():
    def __init__(self, name='NONAME', ambimat='AIR'):
        self.name = name

        self.elements = [] # list of elements comprising this lens [[Element, displacement]]
        self.apertures = {} # dict of apertures ['idx_surface' (aperture comes after this surface), 'z_ap' (absolute position), 'shape', 'radius', 'n_blades']
        self.detector = {} # detector definition ['distance' from last lens, 'sizex', 'sizey', 'pixelpitch', 'pixelsize']

        self.wl = 0.546074 # default wavelength is e-line. unit micrometers
        self.isreverse = False
        self.isbuilt = False

        self.ambimat = ambimat # material between elements and, by default, image and object medium
        self.imagemat = None # optional material behind last lens
        self.objectmat = None # optional material before first lens

        self.idx_image = -1 # index of surface where image is rendered. default -1 (last surface). Used if rendering intermediate images

        # The final data structure
        self.data = {}
        self.num_surfaces = 0

    def add_element(self, ele, dz):
        self.elements.append([ele, dz])

    def build(self, wl=None):
        """This function builds the final data structure"""
        # Data structure definition
        # This dictionary is what the ray tracer actually works thorugh
        self.data['type'] = [] # Surface type for specific functions
        self.data['radius'] = [] # Spherical surface radius
        self.data['asph'] = [] # Apsheric coefficients [[k, a4, a6, ...]]
        self.data['CT_sum'] = [] # Center thickness
        # self.data['lrad'] = [] # Mechanical Lens radius
        self.data['rCA'] = [] # Clear aperture radius
        self.data['decenter'] = [] # Surface lateral position error [[x, y]]
        self.data['tilt'] = [] # Surface tilt error [[X, Y]]
        self.data['n'] = [] # Refractive index
        # self.data['stop'] = [] # Aperture stops [[x,y,z], shape, size]
        # Other parameters
        if wl:
            self.wl = wl

        surfacecount = 0
        # First element of structure is ambient medium plus dummy surface parameters
        if self.objectmat:
            n = glasscatalog.get_n(self.objectmat, self.wl)
        else:
            n = glasscatalog.get_n(self.ambimat, self.wl)
        self.data['type'].append('OBJECT')
        self.data['radius'].append(None)
        self.data['asph'].append(None)
        self.data['CT_sum'].append(None)
        # self.data['lrad'].append(None)
        self.data['rCA'].append(None)
        self.data['decenter'].append(None)
        self.data['tilt'].append(None)
        self.data['n'].append(n)
        surfacecount = surfacecount + 1
        # Add glass surfaces
        CT_sum = 0
        maxrad = 0
        for ele, dz in self.elements:
            num_surfaces = len(ele.data['radius'])
            for s in range(num_surfaces):
                if s == 0:
                    CT_sum = CT_sum + dz
                else:
                    CT_sum = CT_sum + ele.data['CT'][s-1]
                if s == num_surfaces-1:
                    n = glasscatalog.get_n(self.ambimat, self.wl)
                else:
                    n = glasscatalog.get_n(ele.data['material'][s], self.wl)
                self.data['type'].append(ele.data['type'][s])
                self.data['radius'].append(ele.data['radius'][s])
                self.data['asph'].append(ele.data['asph'][s])
                self.data['CT_sum'].append(CT_sum)
                # self.data['lrad'].append(ele.data['lrad'][s])
                self.data['rCA'].append(ele.data['rCA'][s])
                self.data['decenter'].append(ele.data['decenter'][s])
                self.data['tilt'].append(ele.data['tilt'][s])
                self.data['n'].append(n)
                surfacecount = surfacecount + 1
                maxrad = max(maxrad, ele.data['rCA'][s])
        self.isbuilt = True
        self.num_surfaces = surfacecount
        self.maximum_radius = maxrad

        # define detector quad
        z = self.detector['distance'] + self.data['CT_sum'][-1]
        sx = self.detector['sizex']
        sy = self.detector['sizey']
        pp = self.detector['pixelpitch']
        v1 = [-sx*pp/2, -sy*pp/2, z]
        v2 = [sx*pp/2, -sy*pp/2, z]
        v3 = [sx*pp/2, sy*pp/2, z]
        v4 = [-sx*pp/2, sy*pp/2, z]
        self.detector['Quad'] = [v1, v2, v3, v4]

        # convert aperture positions to index
        # TODO: The other way round also?
        for i, ap in self.apertures.items():
            if ap['idx_surface'] is not None:
                # ignore apertures already defined by an index
                continue
            z_ap = ap['z_ap']
            if z_ap <= 0:
                self.apertures[i]['idx_surface'] = 0
                continue
            for j in range(1, self.num_surfaces):
                z_lens = self.data['CT_sum'][j]
                if z_ap > z_lens:
                    self.apertures[i]['idx_surface'] = j
                    

    #def add_aperture(self, aprad, zap, nblades):
    #    self.apertures.append([aprad, zap, nblades])

    #def define_sensor(self, zdet, detsizex, detsizey, relative=True):
    #    dz = 0
    #    if relative:
    #        dz = self.surfaces[-1][2]
    #    self.sensor = [detsizex, detsizey,zdet+dz]

    def change_wavelength(self, wl=0.546074):
        pass

    def reverse_system(self, build=True):
        pass
