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

import os
import re
import numpy as np

from ...raytrace.lenssystem import Lenssystem
from .surfaces import parse_oslo_surface, get_stop
from ..utils import create_elements, get_encoding

# Because of float rounding, values are rounded to these digits for comparison
EPSILON_DIGITS = 3

def _split_opened_oslo(f):
    pass

def split_oslo_file(filename):
    # check encoding of file
    encoding = get_encoding(filename)

    # split he(ader), su(rfaces), fo(oter)
    with open (filename, 'r', encoding=encoding) as f:    
        he, su, fo = _split_opened_oslo(f)

    return he, su, fo

def identify_elements(surf_infos):
    pass

def create_elements(surf_infos, idx_elements, CT_cumulative):
    pass

def load_from_oslo(filename, testmode=False, logfile=None):
    pass