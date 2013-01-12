"""
Tests for splined trajectories.
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

from imusim.trajectories.sampled import SampledTrajectory
from imusim.trajectories.splined import SplinedTrajectory
from imusim.testing.random_data import randomTimeSequence, randomPositionSequence, randomRotationSequence
from imusim.testing.trajectories import checkTrajectory
from imusim.utilities.time_series import TimeSeries

def testSplineTrajectories():
    for i in range(100):
        t = randomTimeSequence()
        p = TimeSeries(t, randomPositionSequence(t))
        r = TimeSeries(t, randomRotationSequence(t))
        T = SplinedTrajectory(SampledTrajectory(p, r), smoothRotations=False)
        yield checkTrajectory, T, p, r
        T = SplinedTrajectory(SampledTrajectory(p, r), smoothRotations=True)
        yield checkTrajectory, T, p, TimeSeries(t, r.values.smoothed())
