"""
Test of Orient-3 IMU model.
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
from imusim.simulation.base import Simulation
from imusim.simulation.calibrators import ScaleAndOffsetCalibrator
from imusim.platforms.imus import Orient3IMU
from imusim.behaviours.imu import BasicIMUBehaviour
from imusim.testing.random_data import RandomTrajectory
from numpy.testing import assert_array_almost_equal
import numpy as np

def testOrient3Model():
    env = Environment()
    calibrator = ScaleAndOffsetCalibrator(env, 1000, 1/100, 10)
    for i in range(3):
        imu = Orient3IMU()
        cal = calibrator.calibrate(imu)
        sim = Simulation(environment=env)
        traj = RandomTrajectory()
        imu.simulation = sim
        imu.trajectory = traj
        sim.time = imu.trajectory.startTime
        BasicIMUBehaviour(imu, 1/256, calibration=cal, initialTime=sim.time)
        sim.run(traj.endTime)
        t = imu.accelerometer.calibratedMeasurements.timestamps
        r = traj.rotation(t)
        a = r.rotateFrame(traj.acceleration(t) -
                env.gravitationalField(traj.position(t), t))
        m = r.rotateFrame(env.magneticField(traj.position(t), t))
        w = r.rotateFrame(traj.rotationalVelocity(t))
        g = env.gravitationalField.nominalMagnitude
        a_within_5g = abs(a) < 5*g
        for i in range(3):
            yield assert_array_almost_equal, \
                np.array([a[i ,a_within_5g[i]]]), \
                np.array([imu.accelerometer.calibratedMeasurements(t)
                    [i, a_within_5g[i]]]), -2
        yield assert_array_almost_equal, m, \
                imu.magnetometer.calibratedMeasurements(t), 0
        yield assert_array_almost_equal, w, \
                imu.gyroscope.calibratedMeasurements(t), 0
