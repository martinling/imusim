"""
Tests for static trajectories
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

from imusim.trajectories.base import StaticTrajectory
from imusim.testing.random_data import randomPosition, randomRotation
from numpy.testing import assert_equal
from imusim.maths.vectors import vector
from imusim.maths.quaternions import QuaternionArray
import numpy as np

def testStaticTrajectory():
    p = randomPosition()
    r = randomRotation()
    s = StaticTrajectory(p, r)
    assert_equal(s.position(0), p)
    assert_equal(s.velocity(0), vector(0,0,0))
    assert_equal(s.acceleration(0), vector(0,0,0))
    assert_equal(s.rotation(0).components, r.components)
    assert_equal(s.rotationalVelocity(0), vector(0,0,0))
    assert_equal(s.rotationalAcceleration(0), vector(0,0,0))
    assert s.startTime == 0
    assert s.endTime == np.inf

    qa = s.rotation((1,2))
    assert type(qa) == QuaternionArray
    assert len(qa) == 2


