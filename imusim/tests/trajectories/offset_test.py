"""
Tests for offset trajectories.
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

from imusim.testing.random_data import randomTimeSequence
from imusim.testing.random_data import randomPosition
from imusim.testing.random_data import randomRotation
from imusim.testing.random_data import RandomTrajectory
from imusim.testing.trajectories import checkTrajectory
from imusim.trajectories.offset import OffsetTrajectory
from imusim.utilities.time_series import TimeSeries

def testOffsetTrajectories():
    for i in range(10):
        T = RandomTrajectory()
        t = T.positionKeyFrames.timestamps
        dt = t[1] - t[0]
        p = T.position(t)
        r = T.rotation(t)
        offset = randomPosition()
        rotationOffset = randomRotation()
        oT = OffsetTrajectory(T, offset, rotationOffset)
        offsetPositions = p + r.rotateVector(offset)
        offsetRotations = r * rotationOffset
        yield checkTrajectory, oT, TimeSeries(t,offsetPositions), TimeSeries(t, offsetRotations)
