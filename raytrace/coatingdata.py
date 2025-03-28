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

"""
In the Element and Lenssystem classes, coating are stored for each surface as a 3-element list. Items are in order:
1) type: Can be 'MIRROR' for ideal mirror, 'DATA' for tabulated AR-coating-data, or None, then Fresnel-reflectance will be computed
2) datasource: filename data source. data must be in column form, column 0 is wavelength in nm, further columns are reflectances in percent
3) dataindex: column index for the data in the file
"""

import numpy as np

class CoatingData():
    def __init__(self, type):
        self.type = type
        self.wl_data = None
        self.reflectance_data = None
        self.filename = None

    def load_coating_csv(self, filename):
        self.filename = filename
        # input must be:
        # column 1: wavelength in nm
        # other columns: refelctance in percent
        data = np.loadtxt(filename, delimiter=',')
        self.wl_data = data[:,0]/1000 # convert to um to follow convention in rest of OptiCore
        self.reflectance_data = data[:, 1:]/100 # convert from percent to abs
        if len(self.reflectance_data.shape) == 1:
            self.reflectance_data = self.reflectance_data[:, np.newaxis]

    def get_R(self, *params):
        # get reflectance
        if self.type == 'DATA':
            return self.get_reflectance_data(params[0], params[1])
        elif self.type == 'MIRROR':
            return 1
        elif self.type == 'FRESNEL_0':
            return self.get_R_fresnel_0(params[0])
        else:
            return self.get_R_fresnel_0(params[0])
        

    def get_reflectance_data(self, wl_um, index):
        # wl outside data range
        if wl_um <= self.wl_data[0]:
            return self.reflectance_data[0, index]
        elif wl_um >= self.wl_data[-1]:
            return self.reflectance_data[-1, index]

        # regular case
        idx_min = np.where(self.wl_data < wl_um)[0][-1]
        idx_max = np.where(self.wl_data > wl_um)[0][0]
        wl_min = self.wl_data[idx_min]
        wl_max = self.wl_data[idx_max]
        r_min = self.reflectance_data[idx_min, index]
        r_max = self.reflectance_data[idx_max, index]

        return r_min + (wl_um - wl_min)/(wl_max - wl_min)*(r_max - r_min)

    def get_R_fresnel_0(self, n_ratio, theta=0):
        # fresnel reflection for unpolarized ligth at 0 degree incidence
        if n_ratio < 1: n_ratio = 1./n_ratio

        R_root = (n_ratio - 1)/(n_ratio + 1)

        return R_root**2

        # TODO: implement numpy-compatible evaluation including angle