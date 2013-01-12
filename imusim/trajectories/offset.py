"""
Trajectories at fixed offsets from other trajectories.
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
from imusim.maths import vectors
from imusim.trajectories.base import PositionTrajectory, RotationTrajectory
from imusim.maths.quaternions import Quaternion
from imusim.utilities.caching import CacheLastValue
from imusim.utilities.documentation import prepend_method_doc

class OffsetTrajectory(PositionTrajectory, RotationTrajectory):
    """
    Trajectory at an offset from another trajectory.

    @ivar parent: the L{AbstractTrajectory} from which this one inherits
    @ivar positionOffset: the offset from the parent trajectory in the parent
        local co-ordinate frame
    @ivar rotationOffset: the rotation offset from the parent trajectory
    """

    def __init__(self, parent, positionOffset=None, rotationOffset=None):
        """
        Initialise trajectory.

        @param parent: the L{AbstractTrajectory} which this trajectory follows.
        @param positionOffset: the offset vector from the parent trajectory in the
            co-ordinate frame of the parent. (3x1 L{np.ndarray})
        @param rotationOffset: the rotational offset from the parent trajectory
            (L{Quaternion})
        """
        self.parent = parent
        if positionOffset is None:
            self.positionOffset = np.zeros((3,1))
        else:
            self.positionOffset = positionOffset
        if rotationOffset is None:
            self.rotationOffset = Quaternion()
        else:
            self.rotationOffset = rotationOffset

    @CacheLastValue()
    def position(self, t):
        p = self.parent.position(t)
        r = self.parent.rotation(t)
        return p + r.rotateVector(self.positionOffset)

    @CacheLastValue()
    def velocity(self, t):
        v = self.parent.velocity(t)
        r = self.parent.rotation(t)
        o = r.rotateVector(self.positionOffset)
        omega = self.parent.rotationalVelocity(t)
        rv = vectors.cross(omega, o)
        return v + rv

    @CacheLastValue()
    @prepend_method_doc(PositionTrajectory)
    def acceleration(self, t):
        """
        Uses the algorithm from Young et al. 2010 'Distributed Estimation
        of Linear Acceleration for Improved Accuracy in Wireless Inertial
        Motion Capture' to calculate linear acceleration from rotational
        velocity and acceleration
        """
        a = self.parent.acceleration(t)
        r = self.parent.rotation(t)
        o = r.rotateVector(self.positionOffset)
        omega = self.parent.rotationalVelocity(t)
        alpha = self.parent.rotationalAcceleration(t)
        lt = vectors.cross(alpha, o)
        lr = vectors.dot(o, omega) * omega - o * vectors.norm(omega)**2
        return a + lt + lr

    @CacheLastValue()
    def rotation(self, t):
        return self.parent.rotation(t) * self.rotationOffset

    def rotationalVelocity(self, t):
        return self.parent.rotationalVelocity(t)

    def rotationalAcceleration(self, t):
        return self.parent.rotationalAcceleration(t)

    @property
    def startTime(self):
        return self.parent.startTime

    @property
    def endTime(self):
        return self.parent.endTime
