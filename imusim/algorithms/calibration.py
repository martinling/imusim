"""
Calibration algorithms for IMU sensors.
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
from scipy.optimize import leastsq
from abc import abstractmethod, ABC
import imusim.maths.vectors as vectors
import numpy as np


class SensorCalibration(ABC):
    """
    Calibration data for a triaxial sensor.
    """
    @abstractmethod
    def apply(measurement):
        """
        Apply calibration to measurement.

        @param measurement: Raw measurements (3x1 L{np.ndarray})
        @return: Calibrated measurements (3x1 L{np.ndarray})
        """
        pass


class ScaleAndOffsetCalibration(SensorCalibration):
    """
    Calibration that corrects constant scale and offset errors.

    Applying this calibration yields:

    M{calibrated = scale * (raw - offset)}

    @ivar scale: Scale factors (3x1 L{np.ndarray})
    @ivar offset: Offset values (3x1 L{np.ndarray})
    """
    def __init__(self, scale, offset):
        """
        Initialise calibration.

        @param scale: Scale factors (3x1 L{np.ndarray})
        @param offset: Offset values (3x1 L{np.ndarray})
        """
        SensorCalibration.__init__(self)

        self.scale = scale
        self.offset = offset

    @staticmethod
    def fromMeasurements(expected, measured):
        """
        Generate calibration by least squares fitting of measured data.

        @param expected: Expected values (3xN L{np.ndarray})
        @param measured: Measured values (3xN L{np.ndarray})

        @return: L{ScaleAndOffsetCalibration} with fitted parameters.

        @raise ValueError: if fit fails to converge.
        """
        params = lambda p: (p[0:3].reshape(3,1), p[3:6].reshape(3,1))
        def error(p):
            cal = ScaleAndOffsetCalibration(*params(p))
            result = cal.apply(measured)
            return vectors.norm(result - expected)
        p0 = np.array([1,1,1,0,0,0])
        p, ier = leastsq(error, p0, ftol=1e-3, maxfev=10000)
        if ier not in [1,2,3,4]:
            raise ValueError("Scale and offset fitting failed.")
        return ScaleAndOffsetCalibration(*params(p))

    def apply(self, measurement):
        return self.scale * (measurement - self.offset)
