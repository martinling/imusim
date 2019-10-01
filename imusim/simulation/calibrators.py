"""
Simulation of calibration procedures for IMUs.
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
from imusim.environment.base import Environment
from imusim.environment.gravity import STANDARD_GRAVITY
from imusim.platforms.sensors import Sensor
from imusim.algorithms.calibration import SensorCalibration, \
        ScaleAndOffsetCalibration
from imusim.trajectories.base import StaticTrajectory, \
        StaticContinuousRotationTrajectory
from imusim.maths.quaternions import Quaternion
import imusim.maths.vectors as vectors
from imusim.simulation.base import Simulation
from imusim.behaviours.imu import BasicIMUBehaviour
from imusim.utilities.documentation import prepend_method_doc
from abc import abstractmethod, ABC
import numpy as np


class Calibrator(ABC):
    """
    A calibration procedure for IMUs.

    @ivar environment: The L{Environment} in which to calibrate.
    """
    def __init__(self, environment):
        """
        Initialise calibrator.

        @param environment: The L{Environment} in which to calibrate.
        """
        self.environment = environment

    @abstractmethod
    def calibrate(imu):
        """
        Calibrate the given L{IMU} device in the simulated environment.

        @return: A dict mapping L{Sensor} objects to L{SensorCalibration}
            objects.

        """
        pass


class ScaleAndOffsetCalibrator(Calibrator):
    """
    Simple calibration procedure that obtains scale and offset factors.

    This procedure simulates testing the device in 18 placements:

        1. Stationary on a flat surface in 6 orthogonal orientations, such
           that gravity is aligned with each direction of each IMU axis.

        2. Stationary on an inclined surface in 6 orthogonal orientations,
           such that magnetic north is aligned with each direction of each
           IMU axis.

        3. On a rotating stage of known angular velocity in 6 orthogonal
           orientations, such that the axis of rotation is aligned with each
           direction of each IMU axis.
    """

    @prepend_method_doc(Calibrator)
    def __init__(self, environment, samples, samplingPeriod,
            rotationalVelocity):
        """
        @param samples: Number of samples to take in each placement (int).
        @param samplingPeriod: Time between samples (float, s).
        @param rotationalVelocity: Angular rate of rotations (float, rad/s).
        """
        Calibrator.__init__(self, environment)
        self.samples = samples
        self.samplingPeriod = samplingPeriod
        self.rotationalVelocity = rotationalVelocity

    def testTrajectories(self):
        position = vectors.vector(0,0,0)
        # Accelerometer test trajectories
        for angles in [
            (0,-90,0),(0,90,0), # X down, X up (pitch -/+90 deg)
            (0,0,90),(0,0,-90), # Y down, Y up (roll +/-90 deg)
            (0,0,0),(0,0,180)]: # Z down, Z up (roll 0/180 deg)
            rotation = Quaternion().setFromEuler(angles)
            yield StaticTrajectory(rotation=rotation)
        # Magnetometer test trajectories
        field = self.environment.magneticField(position, 0)
        N = field / vectors.norm(field)
        E = vectors.vector(-N[1,0], N[0,0], 0)
        E /= vectors.norm(E)
        D = vectors.cross(N, E)
        for vectorset in [
            (N,E,D), (-N,E,-D), # X north, X south
            (E,-N,D), (-E,N,D), # Y north, Y south
            (D,E,-N), (-D,E,N)]: # Z north, Z south
            rotation = Quaternion().setFromVectors(*vectorset)
            yield StaticTrajectory(rotation=rotation)
        # Gyro test trajectories
        for axis in range(3):
            omega = vectors.vector(0,0,0)
            omega[axis] = self.rotationalVelocity
            yield StaticContinuousRotationTrajectory(omega)
            yield StaticContinuousRotationTrajectory(-omega)

    def calibrate(self, imu):
        calibration = {}
        position = vectors.vector(0,0,0)
        expected = {}
        measured = {}
        for sensor in imu.sensors:
            expected[sensor] = []
            measured[sensor] = []
        for trajectory in self.testTrajectories():
            sim = Simulation(environment=self.environment)
            sim.time = 0
            imu.simulation = sim
            imu.trajectory = trajectory
            BasicIMUBehaviour(imu, self.samplingPeriod)
            sim.run(self.samplingPeriod * self.samples,
                    printProgress=False)
            t = sensor.rawMeasurements.timestamps
            for sensor in imu.sensors:
                expected[sensor].append(sensor.trueValues(t))
                measured[sensor].append(sensor.rawMeasurements.values)
        for sensor in imu.sensors:
            calibration[sensor] = ScaleAndOffsetCalibration.fromMeasurements(
                np.hstack(expected[sensor]), np.hstack(measured[sensor]))
        return calibration
