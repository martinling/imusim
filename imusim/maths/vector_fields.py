"""
Classes for modelling 3D vector fields.
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

from abc import ABC, abstractmethod
from scipy import interpolate
from imusim.maths import vectors
from imusim.utilities.documentation import prepend_method_doc
from imusim.maths.natural_neighbour import NaturalNeighbourInterpolatorC
import numpy as np


class VectorField(ABC):
    """
    Base class for vector fields.
    """
    @abstractmethod
    def __call__(self, position, t):
        """
        Evaluate the vector field at the given position(s) and time(s).

        @param position: 3xN L{np.ndarray} of position co-ordinates, in m.
        @param position: Scalar or length N L{np.ndarray} of times in s.

        @return: 3xN L{np.ndarray} of field values.
        """
        pass

    @property
    @abstractmethod
    def nominalValue(self):
        """
        Nominal 3x1 vector value of the field, for use in calibrating sensors.
        """
        pass

    @property
    def nominalMagnitude(self):
        """
        Nominal magnitude of the field, for use in calibrating sensors.
        """
        return vectors.norm(self.nominalValue)


class ConstantVectorField(VectorField):
    """
    A vector field wth a constant value everywhere.
    """

    def __init__(self, value):
        """
        Construct constant field model.

        @param value: Field value as a 3x1 L{np.ndarray}.
        """
        self._value = value

    def __call__(self, position, t):
        result = np.empty_like(position)
        result[:] = self._value
        return result

    @property
    def nominalValue(self):
        return self._value


class InterpolatedVectorField(VectorField):
    """
    A vector field interpolated from sampled values.
    """

    @abstractmethod
    def __init__(self, positions, values):
        """
        Construct vector field interpolator.

        @param positions: 3xN L{np.ndarray} of measurement positions.
        @param values: 3xN L{np.ndarray} of corresponding field values.
        """
        pass


class RBFInterpolatedField(InterpolatedVectorField):
    """
    Field interpolation using radial basis functions.

    Each component of the field is interpolated independently.
    """
    def __init__(self, positions, values):
        x,y,z = positions
        self.components = [interpolate.Rbf(x,y,z,v,function='cubic')
            for v in values]

    def __call__(self, position, t):
        length = np.shape(position)[-1]
        nblocks = min(1, length/500)
        inblocks = np.array_split(position, nblocks, axis=1)
        outblocks = [np.array([np.atleast_1d(c(*(list(ib))))
            for c in self.components]) for ib in inblocks]
        return np.hstack(outblocks)


class NaturalNeighbourInterpolatedField(InterpolatedVectorField):
    """
    Natural Neighbour interpolation of vector fields.

    This is a Python wrapper for the C implementation by Ross Hemsley, described
    in the report "Interpolation on a Magnetic Field", Bristol University, 2009.

    The original code and report are available from:
    http://code.google.com/p/interpolate3d/
    """
    def __init__(self, positions, values):

        valid = vectors.validity(np.vstack((positions,values)))

        x,y,z = positions

        self._baseField = np.median(values[:,valid], axis=1).reshape(3,1)
        deviations = values - self._baseField

        u,v,w = values

        self.imp = NaturalNeighbourInterpolatorC(
            np.ascontiguousarray(x[valid]),
            np.ascontiguousarray(y[valid]),
            np.ascontiguousarray(z[valid]),
            np.ascontiguousarray(u[valid]),
            np.ascontiguousarray(v[valid]),
            np.ascontiguousarray(w[valid]))

    @property
    def nominalValue(self):
        return self._baseField

    def __call__(self, position, t):
        valid = vectors.validity(position)
        result = np.empty_like(position)
        result[:,~valid] = np.nan
        result[:,valid] = self.imp(*(list(position[:,valid])))
        return self._baseField + result
