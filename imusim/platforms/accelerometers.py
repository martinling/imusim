"""
Accelerometer models.
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
from imusim.platforms.sensors import *
from imusim.environment.gravity import STANDARD_GRAVITY
from imusim.utilities.documentation import prepend_method_doc

class Accelerometer(Sensor):
    """
    Base class for accelerometers.
    """

    def trueValues(self, t):
        gravity = self.platform.simulation.environment.gravitationalField
        g = gravity(self.trajectory.position(t), t)
        l = self.trajectory.acceleration(t)
        a = l - g
        return self.trajectory.rotation(t).rotateFrame(a)

class IdealAccelerometer(Accelerometer, IdealSensor):
    """
    An ideal accelerometer.
    """

class IdealGravitySensor(IdealAccelerometer):
    """
    An ideal fictional sensor that measures only acceleration due to gravity.
    """
    def trueValues(self, t):
        gravity = self.platform.simulation.environment.gravitationalField
        g = gravity(self.trajectory.position(t), t)
        return self.trajectory.rotation(t).rotateFrame(-g)

class NoisyAccelerometer(NoisySensor, IdealAccelerometer):
    """
    Accelerometer with additive white Gaussian noise.
    """
    pass

class NoisyTransformedAccelerometer(NoisyTransformedSensor, IdealAccelerometer):
    """
    Accelerometer with affine transform transfer function and Gaussian noise.
    """
    pass

class MMA7260Q(NoisyTransformedSensor, Accelerometer):
    """
    Model of the Freescale MMA7260Q accelerometer.

    @ivar sensitivity: Sensitivity setting, '1.5g', '2g', '4g', or '6g'.
    """

    NOMINAL_SENSITIVITIES = {'1.5g':800E-3, '2g':600E-3,
            '4g':300E-3, '6g':200E-3}
    NOMINAL_OFFSET = 1.65
    MAX_CROSS_AXIS = 0.05

    @prepend_method_doc(Accelerometer)
    def __init__(self, platform, sensitivity, noiseStdDev, rng=None, **kwargs):
        """
        @param sensitivity: Sensitivity setting, '1.5g', '2g', '4g', or '6g'.
        @param noiseStdDev: Standard deviation of measured output voltage noise
            in the system in which the accelerometer is integrated.
        @param rng: L{np.random.RandomState} from which to draw imperfections.
        """
        if rng is None:
            rng = np.random.RandomState()

        sensitivity = MMA7260Q.NOMINAL_SENSITIVITIES[sensitivity] / STANDARD_GRAVITY
        sensitivity *= rng.normal(size=3, loc=1, scale=0.075/3)

        sensitivity = np.diag(sensitivity)
        cross_axis = np.eye(3)
        for s in ((i,j) for i,j in np.ndindex(3, 3) if i != j):
            cross_axis[s] = rng.normal(loc=0, scale=self.MAX_CROSS_AXIS/3)
        transform = np.dot(sensitivity,cross_axis)

        offset = rng.normal(size=(3, 1), loc=self.NOMINAL_OFFSET, scale = 0.165/3)

        NoisyTransformedSensor.__init__(self, platform, noiseStdDev,
                transform, offset, **kwargs)
