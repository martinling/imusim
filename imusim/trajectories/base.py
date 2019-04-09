"""
Base classes for trajectories.
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

from abc import ABC, abstractmethod
from imusim.maths.quaternions import Quaternion, QuaternionArray
from scipy.linalg import expm
import numpy as np


class AbstractTrajectory(ABC):
    """
    Base class of trajectories
    """
    @property
    @abstractmethod
    def startTime(self):
        """
        The first time at which this trajectory is fully defined.
        """
        pass

    @property
    @abstractmethod
    def endTime(self):
        """
        The last time at which this trajectory is fully defined.
        """
        pass


class PositionTrajectory(AbstractTrajectory):
    """Represents a continous trajectory of positions"""

    @abstractmethod
    def position(self, t):
        """
        Return position in the global frame at time t.

        @param t: the time at which to evaluate the position
        @return: the position vector at time t
        """
        pass

    @abstractmethod
    def velocity(self, t):
        """
        Return velocity in the global frame at time t.

        @param t: the time at which to evaluate the velocity
        @return: the velocity vector at time t
        """
        pass

    @abstractmethod
    def acceleration(self, t):
        """
        Return acceleration in the global frame at time t.

        @param t: the time at which to evaluate the acceleration
        @return: the acceleration vector at time t
        """
        pass


class RotationTrajectory(AbstractTrajectory):
    """
    Represents a continuous trajectory of rotations
    """

    @abstractmethod
    def rotation(self,t):
        """
        Return rotation relative to the global frame at time t.

        @param t: the time at which to evaluate the rotation
        @return: the rotation L{Quaternion} at time t
        """
        pass

    @abstractmethod
    def rotationalVelocity(self, t):
        """
        Return rotational velocity relative to the global frame at time t.

        @param t: the time at which to evaluate the rotational velocity
        @return: the rotational velocity vector at time t
        """
        pass

    @abstractmethod
    def rotationalAcceleration(self, t):
        """
        Return rotational acceleration relative to the global frame at time t.

        @param t: the time at which to evaluate the rotational acceleration
        @return: the rotational accelration vector at time t
        """
        pass


class ConstantPositionTrajectory(PositionTrajectory):
    """
    A trajectory with constant position.
    """
    def __init__(self, position=None):
        """
        Initialise the trajectory.

        @param position: the constant position of the object (3x1 L{np.ndarray})
        """
        self.staticPosition = np.zeros((3,1)) if position is None else position

    def position(self, t):
        pos = np.empty((3, len(np.atleast_1d(t))))
        pos[:] = self.staticPosition
        return pos

    def velocity(self,t):
        return np.zeros((3, len(np.atleast_1d(t))))

    def acceleration(self,t):
        return np.zeros((3, len(np.atleast_1d(t))))


class ConstantRotationTrajectory(RotationTrajectory):
    """
    A trajectory with constant rotation.
    """

    def __init__(self, rotation=None):
        """
        Initialise the trajectory.

        @param rotation: the constant rotation of the object (L{Quaternion})
        """
        self.staticRotation = Quaternion() if rotation is None else rotation

    def rotation(self, t):
        if np.shape(t) == ():
            return self.staticRotation
        else:
            rot = np.empty((len(t),4))
            rot[:] = self.staticRotation.components
            return QuaternionArray(rot)

    def rotationalVelocity(self,t):
        return np.zeros((3, len(np.atleast_1d(t))))

    def rotationalAcceleration(self,t):
        return np.zeros((3, len(np.atleast_1d(t))))


class StaticTrajectory(ConstantPositionTrajectory, ConstantRotationTrajectory):
    """
    Represents the trajectory of a static object.
    """
    def __init__(self, position=None, rotation=None):
        """
        Construct static trajectory.

        @param position: the constant position of the object (3x1 L{np.ndarray})
        @param rotation: the constant rotation of the object (L{Quaternion})
        """
        ConstantPositionTrajectory.__init__(self, position)
        ConstantRotationTrajectory.__init__(self, rotation)

    @property
    def startTime(self):
        return 0

    @property
    def endTime(self):
        return np.inf


class ContinuousRotationTrajectory(RotationTrajectory):
    """
    A trajectory with constant rotational velocity.
    """
    def __init__(self, rotationalVelocity, initialRotation=None):
        self.initialRotation = Quaternion() \
            if initialRotation is None else initialRotation
        self._rotationalVelocity = rotationalVelocity
        w = rotationalVelocity
        W = np.asmatrix([[ 0,      -w[2,0],  w[1,0]],
                         [ w[2,0],  0,      -w[0,0]],
                         [-w[1,0],  w[0,0],  0     ]])
        self.qW = Quaternion()
        self.qW.setFromMatrix(expm(W))

    def rotation(self, t):
        if np.ndim(t) == 0:
            return self.qW ** t
        else:
            return QuaternionArray([self.qW ** ti for ti in t])

    def rotationalVelocity(self, t):
        omega = np.empty((3, len(np.atleast_1d(t))))
        omega[:] = self._rotationalVelocity
        return omega

    def rotationalAcceleration(self, t):
        return np.zeros((3, len(np.atleast_1d(t))))


class StaticContinuousRotationTrajectory(ContinuousRotationTrajectory,
        ConstantPositionTrajectory):
    """
    A trajectory with constant position and rotational velocity.
    """
    def __init__(self, rotationalVelocity, initialRotation=None,
            position=None):
        ConstantPositionTrajectory.__init__(self, position)
        ContinuousRotationTrajectory.__init__(self, rotationalVelocity,
                initialRotation)

    @property
    def startTime(self):
        return 0

    @property
    def endTime(self):
        return np.inf
