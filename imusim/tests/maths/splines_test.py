"""
Tests for scalar splines.
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

from imusim.maths.splines import Spline, UnivariateSpline, PartialInputSpline
from imusim.testing.random_data import randomValidity, invalidate
from nose.tools import raises
from numpy.testing import assert_almost_equal
from numpy import random
import numpy as np

@raises(Spline.InsufficientPointsError)
def testInsufficientPointsError():
    x = [1,2,3]
    y = x
    spline = UnivariateSpline(x, y, order=5)

def testInterpolatingSpline():
    x = np.linspace(0,10,100)
    y = x**3
    yp = 3*x**2
    ypp = 6*x

    spline = UnivariateSpline(x,y, order=3)
    s = (x>spline.validFrom) & (x < spline.validTo)
    assert_almost_equal(spline(x)[s], y[s])
    assert_almost_equal(spline(x, 1)[s], yp[s])
    assert_almost_equal(spline(x, 2)[s], ypp[s])

def testPartialInterpolatingSpline():
    x = np.linspace(0,10,1000)
    x = invalidate(x, randomValidity(x))
    y = x**3
    yp = 3*x**2
    ypp = 6*x

    spline = PartialInputSpline(x,y, order=3)
    valid = np.array(
        np.sum(((x >= spline.validFrom) & (x <= spline.validTo) for spline in
            spline.splines)), dtype=bool)
    assert_almost_equal(spline(x)[valid], y[valid])
    assert_almost_equal(spline(x, 1)[valid], yp[valid])
    assert_almost_equal(spline(x, 2)[valid], ypp[valid])

def testSmootingSpline():
    x = np.linspace(-10,10,1000)
    y = x**3
    random.seed(0)

    spline = UnivariateSpline(x, y + random.normal(scale=0.1, size=1000),
            order=3, stddev=0.1)

    s = (x>spline.validFrom) & (x < spline.validTo)
    assert_almost_equal(spline(x)[s], y[s], 1)



