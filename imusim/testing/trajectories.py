"""
Utilities for testing trajectories.
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
from imusim.testing.quaternions import assertQuaternionAlmostEqual
from imusim.maths.quaternions import QuaternionArray
from imusim.testing.vectors import assert_vectors_correlated
from imusim.utilities.time_series import TimeSeries
import numpy as np

def checkTrajectory(T, truePositions, trueRotations):
    """
    Check the outputs of a trajectory model agree with truth values.

    @param T: Trajectory to check.
    @param truePositions: L{TimeSeries} of true position values.
    @param trueRotations: L{TimeSeries} of true rotation values.
    """

    # Get time indices at which position comparisons valid
    t = truePositions.timestamps
    validity = (t >= T.startTime) & (t <= T.endTime)
    t = t[validity]
    dt = np.gradient(t)
    p = truePositions.values[:,validity]

    # Check position
    assert_vectors_correlated(T.position(t), p)

    # Check velocity
    v = np.array(list(map(np.gradient, p))) / dt
    assert_vectors_correlated(T.velocity(t[2:-2]), v[:,2:-2])

    # Check acceleration
    a = np.array(list(map(np.gradient, v))) / dt
    assert_vectors_correlated(T.acceleration(t[4:-4]), a[:,4:-4])

    # Get time indices at which rotation comparisons valid
    t = trueRotations.timestamps
    validity = (t >= T.startTime) & (t <= T.endTime)
    t = t[validity]
    r = trueRotations.values[validity]

    # Check rotation
    assertQuaternionAlmostEqual(T.rotation(t), r, tol=0.05)

    # Check angular velocity
    r, lastR = r[1:], r[:-1]
    t, dt = t[1:], np.diff(t)
    diffOmega = (2 * (r - lastR) * lastR.conjugate).array.T[1:] / dt
    trajOmega = T.rotationalVelocity(t - dt/2)
    assert_vectors_correlated(trajOmega[:,2:-2], diffOmega[:,2:-2])

    # Check angular acceleration
    diffAlpha = np.array(list(map(np.gradient, diffOmega))) / dt
    trajAlpha = T.rotationalAcceleration(t - dt/2)
    assert_vectors_correlated(trajAlpha[:,4:-4], diffAlpha[:,4:-4])
