"""
Tests for column vector utilities.
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
from numpy.testing import assert_almost_equal
from imusim.maths import vectors

e1 = np.array([[1.0,0,0]]).T
e2 = np.array([[0,1.0,0]]).T
e3 = np.array([[0,0,1.0]]).T

basisVectors = (e1,e2,e3)

def testNorm():
    for v in basisVectors:
        assert_almost_equal(vectors.norm(v), 1)

def testCross():
    assert_almost_equal(vectors.cross(e1, e2), e3)
    assert_almost_equal(vectors.cross(e1, e3), -e2)
    assert_almost_equal(vectors.cross(e2, e3), e1)

def testDot():
    for v1 in basisVectors:
        for v2 in basisVectors:
            if v1 is v2:
                assert_almost_equal(vectors.dot(v1, v2), 1)
            else:
                assert_almost_equal(vectors.dot(v1, v2), 0)
