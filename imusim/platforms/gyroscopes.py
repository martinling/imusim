"""
Gyroscope models.
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
import math
from imusim.platforms.sensors import Sensor, IdealSensor, NoisySensor, \
    NoisyTransformedSensor
from imusim.platforms.accelerometers import Accelerometer
from imusim.environment.gravity import STANDARD_GRAVITY
from imusim.utilities.documentation import prepend_method_doc

class Gyroscope(Sensor):
    """
    Base class for gyroscopes.
    """

    def trueValues(self, t):
        omega = self.trajectory.rotationalVelocity(t)
        return self.trajectory.rotation(t).rotateFrame(omega)

class IdealGyroscope(Gyroscope, IdealSensor):
    """
    An ideal gyroscope.
    """
    pass

class NoisyGyroscope(NoisySensor, IdealGyroscope):
    """
    Gyroscope with additive white Gaussian noise.
    """
    pass

class NoisyTransformedGyroscope(NoisyTransformedSensor, IdealGyroscope):
    """
    Gyroscope with affine transform transfer function and Gaussian noise.
    """
    pass

class ADXRS300(NoisyTransformedSensor, Gyroscope):
    """
    Model of the Analog devices ADXRS300 gyroscope with enhanced range.

    For range extension method see Analog Devices application note AN625.

    The model assumes that the output of the gyroscope is converted to a
    3.3V voltage domain using a 1% tolerance resistor divider.
    """

    NOMINAL_SENSITIVITY = {'300deg/s':5E-3, '600deg/s':2.5E-3,
            '900deg/s':1.67E-3, '1200deg/s':1.25E-3} # mV/deg/s
    NOMINAL_OFFSET = 2.5 # V
    CROSS_AXIS = 0.05 # 5% cross axis
    NOMINAL_ACCELERATION_SENSITIVITY = np.radians(0.2) # rad/s/g

    @prepend_method_doc(Gyroscope)
    def __init__(self, platform, noiseStdDev, sensitivity='300deg/s',
            rng=None, **kwargs):
        """
        @param sensitivity: Sensitivity, '300deg/s', '600deg/s', '900deg/s'
            or '1200deg/s'.
        @param noiseStdDev: Standard deviation of measured output voltage noise
            in the system in which the gyroscope is integrated.
        @param rng: L{np.random.RandomState} from which to draw imperfections.
        """
        if rng is None:
            rng = np.random.RandomState()

        sensitivity = self.NOMINAL_SENSITIVITY[sensitivity] * 180 / np.pi
        sensitivity *= rng.normal(size=3, loc=1, scale=0.08/3)

        sensitivity = np.diag(sensitivity)
        cross_axis = np.eye(3)
        for s in ((i,j) for i,j in np.ndindex(3, 3) if i != j):
            cross_axis[s] = rng.normal(loc=0, scale=self.CROSS_AXIS/3)
        transform = np.dot(sensitivity, cross_axis)

        offset = rng.normal(size=(3, 1), loc=self.NOMINAL_OFFSET, scale=0.02/3)

        # Conversion from 5V domain to 3.3V
        # Errors in 1% divider resistors have negligible effect compared to
        # 8% error in gyroscope sensitivity
        transform *= (3.3/5)
        offset *= (3.3/5)

        NoisyTransformedSensor.__init__(self, platform, noiseStdDev,
                transform, offset, **kwargs)

        self._accelSensitivity = rng.normal(size=(3, 1), loc=0,
                scale=self.NOMINAL_ACCELERATION_SENSITIVITY/3)

    def sensedVoltages(self, t):
        # Measure acceleration to model acceleration sensitivity
        gravity = self.platform.simulation.environment.gravitationalField
        g = gravity(self.trajectory.position(t), t)
        l = self.trajectory.acceleration(t)
        a = l - g
        aLocal = self.trajectory.rotation(t).rotateFrame(a)
        return self._transform.apply(self.trueValues(t) \
                + aLocal * self._accelSensitivity)
