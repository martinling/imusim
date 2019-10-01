"""
IMU hardware platform models.
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
from imusim.platforms.base import Platform
from imusim.platforms.sensors import Sensor
from imusim.platforms.accelerometers import Accelerometer, \
        IdealAccelerometer, IdealGravitySensor, MMA7260Q
from imusim.platforms.magnetometers import Magnetometer, \
        IdealMagnetometer, HMC105x
from imusim.platforms.gyroscopes import Gyroscope, \
        IdealGyroscope, ADXRS300
from imusim.platforms.adcs import ADC, IdealADC, QuantisingADC
from imusim.platforms.timers import Timer, IdealTimer, ParametricTimer
from imusim.platforms.radios import Radio, IdealRadio
from imusim.maths.vectors import vector
from imusim.utilities.documentation import prepend_method_doc


class IMU(Platform):
    """
    An IMU hardware platform with one or more sensors.
    """
    @property
    @abstractmethod
    def sensors(self):
        """
        List of L{Sensor} objects on the platform.
        """
        pass


class StandardIMU(IMU):
    """
    An IMU with accelerometer, magnetometer, gyroscope, ADC, Timer and radio.

    More complex designs may not fit this structure, but it provides a useful
    interface which simplifies repeating experiments across multiple IMU types.

    @ivar accelerometer: The L{Accelerometer} on the platform.
    @ivar magnetometer: The L{Magnetometer} on the platform.
    @ivar gyroscope: The L{Gyroscope} on the platform.
    @ivar adc: The L{ADC} used to sample the sensors.
    @ivar timer: The L{Timer} used for timekeeping on the platform.
    @ivar radio: The L{Radio} on the platform.
    """

    @property
    def sensors(self):
        return [self.accelerometer, self.magnetometer, self.gyroscope]

    @property
    def components(self):
        return self.sensors + [self.adc, self.timer, self.radio]


class IdealIMU(StandardIMU):
    """
    An IMU with idealised models for all components.
    """

    def __init__(self, simulation=None, trajectory=None):
        self.accelerometer = IdealAccelerometer(self)
        self.magnetometer = IdealMagnetometer(self)
        self.gyroscope = IdealGyroscope(self)
        self.adc = IdealADC(self)
        self.timer = IdealTimer(self)
        self.radio = IdealRadio(self)
        StandardIMU.__init__(self, simulation, trajectory)


class MagicIMU(StandardIMU):
    """
    An IMU with idealised components including a fictional gravity sensor.
    """
    def __init__(self, simulation=None, trajectory=None):
        self.accelerometer = IdealGravitySensor(self)
        self.magnetometer = IdealMagnetometer(self)
        self.gyroscope = IdealGyroscope(self)
        self.adc = IdealADC(self)
        self.timer = IdealTimer(self)
        self.radio = IdealRadio(self)
        StandardIMU.__init__(self, simulation, trajectory)


class Orient3IMU(StandardIMU):
    """
    Simple model of the Orient-3 IMU from the University of Edinburgh.

    The Orient-3 is an updated version of the Orient-2 design described in
    A. Young, M. Ling and D. K. Arvind, "Orient-2: A Realtime Wireless Posture
    Tracking System Using Local Orientation Estimation", in Proc. 4th Workshop
    on Embedded Network Sensors, pp 53-57, ACM, 2007.

    This model still excludes a number of factors such as sampling delays, and
    currently includes an ideal radio model, but provides reasonably realistic
    sensor data. The sensor models are based on a combination of empirical
    characterisation and datasheet parameters.
    """

    @prepend_method_doc(Platform)
    def __init__(self, simulation=None, trajectory=None, rng=None):
        """
        @param rng: L{np.random.RandomState} from which to draw imperfections.
        """
        self.adc = QuantisingADC(self, bits=12, vref=1.65)
        adc_lsb = 3.3/2**self.adc._bits
        self.timer = ParametricTimer(self, frequency=32768, freqError=30, rng=rng)
        self.accelerometer = MMA7260Q(self, sensitivity='6g',
                noiseStdDev=15.3*adc_lsb, rng=rng)
        self.magnetometer = HMC105x(self, noiseStdDev=8.1*adc_lsb, rng=rng)
        self.gyroscope = ADXRS300(self, sensitivity='1200deg/s',
                noiseStdDev=1*adc_lsb, rng=rng)
        self.radio = IdealRadio(self)
        StandardIMU.__init__(self, simulation, trajectory)
