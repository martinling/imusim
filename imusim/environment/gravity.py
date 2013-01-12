"""
Gravitational field models.
"""
# Copyright (C) 2009-2011 University of Edinburgh
#
# This file is part of IMUSim.
#
# IMUSim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# IMUSim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with IMUSim.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from imusim.maths.vectors import vector
from imusim.maths.vector_fields import ConstantVectorField

""" Standard Earth gravity in m/s^2. """
STANDARD_GRAVITY = 9.81
"""Standard gravitational field magnitdue in M{m/s^2}"""

class ConstantGravitationalField(ConstantVectorField):
    """
    Constant gravitational field model.

    The field is modelled as a parallel vector field with a constant
    magnitude. The field is non-zero only in the z (Down) axis.
    """

    def __init__(self, magnitude=STANDARD_GRAVITY):
        """
        Construct field model.

        @param magnitude: Magnitude of field in m/s^2. Default is standard
            Earth gravity.
        """
        ConstantVectorField.__init__(self, vector(0, 0, magnitude))

    @property
    def nominalValue(self):
        return self._value.copy()

    def __call__(self, position, t):
        field = np.empty_like(position)
        field[:] = self._value
        return field
