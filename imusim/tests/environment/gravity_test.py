"""
Tests for gravitational field environment.
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

from imusim.environment.gravity import ConstantGravitationalField
import numpy as np
from numpy.testing import assert_equal

def testGravitationalField():
    f = ConstantGravitationalField()

    assert_equal(f.nominalMagnitude,9.81)
    assert_equal(f(np.zeros((3,1)), 0),f.nominalValue)
