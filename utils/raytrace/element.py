class Element():
    def __init__(self, name='NONAME'):
        self.name = name
        self.data = {}
        self.data['type'] = [] # Surface type for specific functions
        self.data['radius'] = [] # Spherical surface radius
        self.data['asph'] = [] # Apsheric coefficients [[k, a4, a6, ...]]
        self.data['CT'] = [] # Center thickness
        # self.data['lrad'] = [] # Mechanical Lens radius
        self.data['rCA'] = [] # Clear aperture radius
        self.data['decenter'] = [] # Surface lateral position error [[x, y]]
        self.data['tilt'] = [] # Surface tilt error [[X, Y]]
        self.data['material'] = [] # Glass type specifiers

        self.isreverse = False
        self.isdummy = False # Dummy lenses for placing e.g. an aperture

    def reverse(self):
        for key in self.data:
            if key == 'radius':
                self.data[key] = [-1.*r for r in self.data[key][::-1]]
            elif key == 'asph':
                self.data[key] = [[e[0]] + [-1.*e[i] for i in range(1,len(e))] for e in self.data[key][::-1]] # conical constant remains unchanged
            elif key == 'tilt':
                self.data[key] = [[-1.*e[0], -1.*e[1]] for e in self.data[key][::-1]]
            else:
                self.data[key] = self.data[key][::-1]
        self.isreverse =  not self.isreverse