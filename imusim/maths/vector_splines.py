"""
Spline fitting of vector data.
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

from __future__ import division
from .splines import Spline, UnivariateSpline, PartialInputSpline
import numpy as np
from imusim.maths import vectors

class UnivariateVectorSpline(Spline):
    """
    Model of a vector function of a single variable using spline fitting.

    The implementation uses an independent L{UnivariateSpline} for each
    component of the vector function.
    """
    def __init__(self, x, y, **kwargs):
        """
        Construct vector spline.

        @param x: Length N, monotonically increasing sequence of input values.
        @param y: MxN L{np.ndarray} of output values.
        @param kwargs: Additional parameters passed to L{UnivariateSpline}.
        """
        self.splines = [UnivariateSpline(x, yi, **kwargs) for yi in y]

    def __call__(self,x,n=0):
        """
        Evaluate the n-th derivative of the function at input value(s) x.
        """
        return np.vstack([s(x, n) for s in self.splines])

    @property
    def validFrom(self):
        return max([s.validFrom for s in self.splines])

    @property
    def validTo(self):
        return min([s.validTo for s in self.splines])

class PartialInputVectorSpline(PartialInputSpline):
    """
    Piecewise vector spline allowing for undefined regions in output domain.
    """
    _splineClass = UnivariateVectorSpline

    def __init__(self, x, y, **kwargs):
        PartialInputSpline.__init__(self, x, y, **kwargs)
        self._dims = np.shape(y)[0]

    def _validity(self, y):
        return vectors.validity(y)

    def _output(self, x, conditions, results, undefined):
        out = np.empty((self._dims, len(np.atleast_1d(x))))
        for cond, result in zip(conditions, results):
            out[:,cond] = result
        out[:,undefined] = np.nan
        return out
