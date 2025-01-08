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

def get_N1_sqsurface(N2, dshape):
    # passing N2 and calculating N1, not the other way, for historical reasons
    # TODO: Might be a good idea to restrict N2 to odd numbers? Particularly for dshape cases to retain square segments. Consider depending on results
    if dshape:
        N1 = N2//2 + 1
    else:
        N1 = N2
    return N1