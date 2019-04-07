"""
Behaviours for IMU devices.
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
from abc import ABC
from imusim.utilities.time_series import TimeSeries
from imusim.behaviours.timing import TimerMultiplexer, VirtualTimer
from imusim.behaviours.sampling import PeriodicSampler
from imusim.algorithms.calibration import SensorCalibration
from imusim.platforms.sensors import Sensor
from imusim.algorithms.orientation import OrientationFilter
from imusim.platforms.imus import IMU


class BasicIMUBehaviour(ABC):
    """
    Basic behaviour for an IMU that performs periodic sampling.

    @ivar imu: L{IMU} on which the behaviour executes.
    @ivar samplingPeriod: Interval at which to sample and process (float).
    @ivar calibration: Mapping from IMU L{Sensor}s to L{SensorCalibration}s
        (dict).
    @ivar filter: L{OrientationFilter} to update with each set of samples.
    @ivar sampleCallback: Function to call after each sample. The behaviour
        object will be passed as a single argument.
    @ivar timerMux: L{TimerMultiplexer} for IMU timer.
    """
    def __init__(self, imu, samplingPeriod, calibration=None, filter=None,
            sampleCallback=None, initialTime=0):
        """
        Initialise IMU behaviour.

        @param imu: L{IMU} on which the behaviour executes.
        @param samplingPeriod: Interval at which to sample and process (float).
        @param calibration: Mapping from IMU L{Sensor}s to
            L{SensorCalibration}s (dict).
        @param filter: L{OrientationFilter} to update with each set of samples.
        @param sampleCallback: Function to call after each sample. The
            behaviour object will be passed as a single argument.
        @param initialTime: Initial time for local timekeeping.
        """

        self.imu = imu
        self.samplingPeriod = samplingPeriod
        self.calibration = calibration
        self.filter = filter

        for sensor in imu.sensors:
            sensor.rawMeasurements = TimeSeries()
            if self.calibration is not None:
                sensor.calibratedMeasurements = TimeSeries()

        self.sampleCallback = sampleCallback

        self.timerMux = TimerMultiplexer(imu.timer)

        PeriodicSampler(imu.adc, imu.sensors, VirtualTimer(self.timerMux),
                samplingPeriod, self._handleSample)

        self._time = initialTime

    def _handleSample(self, *measurements):
        measurements = list(measurements)
        self._time += self.samplingPeriod
        for sensor, measurement in zip(self.imu.sensors, measurements):
            sensor.rawMeasurements.add(self._time, measurement)
            if self.calibration is not None:
                measurement = self.calibration[sensor].apply(measurement)
                sensor.calibratedMeasurements.add(self._time, measurement)
        if self.filter is not None:
            self.filter(*measurements + [self._time])
        if self.sampleCallback != None:
            self.sampleCallback(self)
