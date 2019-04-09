"""
Base classes for sensor models.
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

from abc import abstractmethod
from imusim.platforms.base import Component
from imusim.maths.transforms import AffineTransform
from imusim.trajectories.offset import OffsetTrajectory
from imusim.simulation.base import Simulation
from imusim.utilities.documentation import prepend_method_doc
from imusim.maths.vectors import vector
from imusim.maths.quaternions import Quaternion
import numpy as np


class Sensor(Component):
    """
    Base class for all IMU sensor classes.

    @ivar platform: L{Platform} this sensor is attached to.
    @ivar trajectory: L{OffsetTrajectory} followed by this sensor.
    """
    def __init__(self, platform, positionOffset=vector(0,0,0),
            rotationOffset=Quaternion(1,0,0,0)):
        """
        Initialise sensor.

        @param platform: L{Platform} this sensor will be attached to.
        @param positionOffset: Position offset vector of this sensor
            in the local co-ordinate frame of its host platform.
        @param rotationOffset: Rotation offset quaternion of this sensor
            in the local co-ordinate frame of its host platform.
        """
        self._positionOffset = positionOffset
        self._rotationOffset = rotationOffset
        Component.__init__(self, platform)

    def _trajectoryChange(self):
        parentTrajectory = self.platform.trajectory
        if parentTrajectory is None:
            self.trajectory = None
        else:
            self.trajectory = OffsetTrajectory(parentTrajectory,
                self._positionOffset, self._rotationOffset)

    @abstractmethod
    def trueValues(self, t):
        """
        Get true values of the intended sensed quantity at given times.

        @param t: time(s) at which to measure.

        @return: 3xN L{np.ndarray} of true values, where N = len(t).
        """
        pass

    @abstractmethod
    def sensedVoltages(self, t):
        """
        Get output voltages at given times.

        This function should include all deterministic effects in the sensor,
        but not random noise.

        @param t: time(s) at which to measure.

        @return: 3xN L{np.ndarray} of output voltages, where N = len(t).
        """
        pass

    @abstractmethod
    def noiseVoltages(self, t):
        """
        Generate non-deterministic noise to be added to output voltages.

        @param t: time(s) at which measurements were taken, len(t) = N.

        @return: 3xN L{np.ndarray} of noise voltages.
        """
        pass

    def voltages(self, t):
        """
        Generate output voltages of the sensor at given times, including noise.
        """
        return self.sensedVoltages(t) + self.noiseVoltages(t)


class IdealSensor(Sensor):
    """
    An ideal sensor with unity transfer function and no noise.
    """

    def sensedVoltages(self, t):
        return self.trueValues(t)

    def noiseVoltages(self, t):
        return vector(0,0,0)


class NoisySensor(Sensor):
    """
    An sensor with additive white Gaussian noise.
    """

    @prepend_method_doc(Sensor)
    def __init__(self, platform, noiseStdDev, **kwargs):
        """
        @param noiseStdDev: Standard deviation of measurement noise.
        """
        Sensor.__init__(self, platform, **kwargs)
        self._noiseStdDev = noiseStdDev

    def _simulationChange(self):
        Sensor._simulationChange(self)
        if self.platform.simulation is None:
            self._rng = None
        else:
            self._rng = self.platform.simulation.subrng()

    def noiseVoltages(self, t):
        return self._rng.normal(size=(3,1), scale=self._noiseStdDev)


class TransformedSensor(Sensor):
    """
    An sensor with an affine transform transfer function.
    """

    @prepend_method_doc(Sensor)
    def __init__(self, platform, matrix=np.eye(3), offset=np.zeros((3,1)),
            **kwargs):
        """
        @param matrix: 3x3 matrix to transform true values.
        @param offset: 3x1 offset to be added after matrix multiplication.
        """
        Sensor.__init__(self, platform, **kwargs)
        self._transform = AffineTransform(transform=matrix, translation=offset)

    def sensedVoltages(self, t):
        return self._transform.apply(self.trueValues(t))


class NoisyTransformedSensor(NoisySensor, TransformedSensor):
    """
    A sensor with an affine transform transfer function and Gaussian noise.
    """
    @prepend_method_doc(NoisySensor)
    def __init__(self, platform, noiseStdDev, transform, offset, **kwargs):
        """
        @param transform: 3x3 transform matrix to transform true values.
        @param offset: 3x1 offset to be added after matrix multiplication.
        """
        NoisySensor.__init__(self, platform, noiseStdDev, **kwargs)
        TransformedSensor.__init__(self, platform, transform, offset, **kwargs)
