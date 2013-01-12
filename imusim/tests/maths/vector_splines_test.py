"""
Tests for vector splines.
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

from imusim.maths.splines import Spline
from imusim.maths.vector_splines import UnivariateVectorSpline
from imusim.maths.vector_splines import PartialInputVectorSpline
from imusim.testing.random_data import randomTimeSequence
from imusim.testing.random_data import randomPositionSequence
from imusim.testing.random_data import randomValidity, invalidate
from numpy.testing import assert_almost_equal
import numpy as np

def checkVectorSpline(splineClass, x, y, order):
    try:
        spline = splineClass(x, y, order=order)
        ys = spline(x)
        splines = spline.splines if isinstance(spline, PartialInputVectorSpline) else [spline]
        valid = np.array(
            np.sum(((x >= spline.validFrom) & (x <= spline.validTo) for spline in splines)),
            dtype=bool)
        assert_almost_equal(y[:,valid], ys[:,valid])
    except Spline.InsufficientPointsError:
        pass

def testVectorSpline():
    for order in (1, 3, 5):
        for i in range(10):
            x = randomTimeSequence()
            y = randomPositionSequence(x)
            yield checkVectorSpline, UnivariateVectorSpline, x, y, order
            validity = randomValidity(x)
            yn = invalidate(y, validity)
            yield checkVectorSpline, PartialInputVectorSpline, x, yn, order
