"""
Magnetometer models.
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
from imusim.utilities.documentation import prepend_method_doc

class Magnetometer(Sensor):
    """
    Base class for magnetometers.
    """

    def trueValues(self, t):
        b = self.platform.simulation.environment.magneticField(
                self.trajectory.position(t), t)
        return self.trajectory.rotation(t).rotateFrame(b)

class IdealMagnetometer(Magnetometer, IdealSensor):
    """
    An ideal magnetometer.
    """
    pass

class NoisyMagnetometer(NoisySensor, IdealMagnetometer):
    """
    Magnetometer with additive white Gaussian noise.
    """
    pass

class NoisyTransformedMagnetometer(NoisyTransformedSensor, IdealMagnetometer):
    """
    Magnetometer with an affine transform transfer function and Gaussian noise.
    """
    pass

class HMC105x(NoisyTransformedSensor, Magnetometer):
    """
    Model of the Honeywell HMC105x magnetometer as used in the Orient-3 IMU.

    B{N.B. This model includes an INA321 instrument amplifier circuit.}

    The HMC105x has a sensitivity of 1mV/V/Gauss +- 20%. IMUSim uses Teslas for
    magnetic field measurements, where 1 Tesla = 10,000 Gauss.

    The Orient-3 design uses an INA321 precision instrument amplifier with
    a gain of 380+-1% to amplify the magnetometer signal prior to sampling.
    """

    VOLTS_PER_VOLT_PER_GAUSS = 1E-3
    GAUSS_PER_TESLA = 10000
    VDD = 3.3
    SENSITIVITY = VOLTS_PER_VOLT_PER_GAUSS * GAUSS_PER_TESLA * VDD
    MAX_SENSITIVITY_ERROR = 0.2 # 20%
    INA321_NOMINAL_GAIN = 380
    BRIDGE_OFFSET = 0.5E-3 * VDD
    INA321_INPUT_OFFSET = 200E-6
    CROSS_AXIS = 0.05

    @prepend_method_doc(Magnetometer)
    def __init__(self, platform, noiseStdDev, rng=None, **kwargs):
        """
        @param noiseStdDev: Standard deviation of measured output voltage
            noise in the system in which the magnetometer is integrated.
        @param rng: L{np.random.RandomState} from which to draw imperfections.
        """

        if rng is None:
            rng = np.random.RandomState()

        sensitivity = self.SENSITIVITY
        sensitivity *= self.INA321_NOMINAL_GAIN
        sensitivity *= rng.normal(size=3, loc=1,
                scale=self.MAX_SENSITIVITY_ERROR/3)

        sensitivity = np.diag(sensitivity)
        cross_axis = np.eye(3)
        for s in ((i,j) for i,j in np.ndindex(3, 3) if i != j):
            cross_axis[s] = rng.normal(loc=0, scale=self.CROSS_AXIS/3)
        transform = np.dot(sensitivity,cross_axis)

        offset = self.VDD / 2
        offset += rng.normal(size=(3, 1), loc=0,
                scale=self.INA321_INPUT_OFFSET * self.INA321_NOMINAL_GAIN/3)
        offset += rng.normal(size=(3, 1), loc=0,
                scale=self.BRIDGE_OFFSET * self.INA321_NOMINAL_GAIN/3)

        NoisyTransformedSensor.__init__(self, platform, noiseStdDev,
                transform, offset, **kwargs)
