"""
Tests for orientation tracking algorithms.
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
from imusim.testing.random_data import RandomRotationTrajectory
from imusim.testing.quaternions import assert_quaternions_correlated, assertMaximumErrorAngle
from imusim.testing.inspection import getImplementations
from imusim.algorithms import orientation
from imusim.maths.quaternions import Quaternion, QuaternionArray
from imusim.trajectories.rigid_body import Joint
from imusim.maths.vectors import vector

dt = 0.01
GRAVITY = vector(0,0,1)
NORTH = vector(1,0,0)

filterParameters = {
    orientation.OrientCF : {
        'k' : 1,
        'aT' : 1},
    orientation.BachmannCF : {
        'k' : 1},
    orientation.DistLinAccelCF : {
        'k' : 1,
        'joint' : Joint(None),
        'offset' : np.zeros((3,1)),
        'gravFieldReference' : GRAVITY,
        'magFieldReference' : NORTH
        },
    orientation.YunEKF : {
        'initialCovariance' : np.asmatrix(
            np.diag([0.01]*3 + [0.0001]*4)),
        'measurementCovariance' : np.asmatrix(
            np.diag([0.01]*3 + [0.0001]*4)),
        'D' : 50,
        'tau' : 0.5},
    }

def runOrientationFilter(filter, trajectory):
    t = np.arange(trajectory.startTime + dt, trajectory.endTime, dt)
    trueOrientations = trajectory.rotation(t)
    gyroSamples = trueOrientations.rotateFrame(
            trajectory.rotationalVelocity(t))
    accelSamples = trueOrientations.rotateFrame(-GRAVITY)
    magSamples = trueOrientations.rotateFrame(NORTH)
    estimatedOrientations = QuaternionArray(np.empty_like(
        trueOrientations.array))

    for i in range(gyroSamples.shape[1]):
        accel = accelSamples[:,i:i+1]
        mag = magSamples[:,i:i+1]
        gyro = gyroSamples[:,i:i+1]
        estimatedOrientations[i] = filter.rotation.latestValue
        filter(accel, mag, gyro, t[i])

    return trueOrientations, estimatedOrientations

def checkOrientationFilter(filter, trajectory):
    trueOrientations, estimatedOrientations = runOrientationFilter(filter,
            trajectory)
    assert_quaternions_correlated(estimatedOrientations, trueOrientations)

def checkDynamicConvergence(filter, trajectory):
    trueOrientations, estimatedOrientations =  runOrientationFilter(filter,
            trajectory)
    assert_quaternions_correlated(trueOrientations[-100:],
            estimatedOrientations[-100:])

def checkStaticConvergence(filterClass, performTest=True):
    """
    Test that an orientation estimation algorithm converges within a reasonable
    number of samples
    """
    gyro = np.zeros((3,1))
    accel = -GRAVITY
    mag = NORTH
    SAMPLES = 500

    initialRotation = Quaternion.fromEuler((60,45,0))

    filter = filterClass(initialRotation = initialRotation,
            initialTime=0,
            initialRotationalVelocity = np.zeros((3,1)),
            **filterParameters.get(filterClass, {}))

    estimatedQuaternions = QuaternionArray(np.empty((SAMPLES,4)))

    for i in range(SAMPLES):
        filter(accel, mag, gyro, dt*(1+i))
        estimatedQuaternions[i] = filter.rotation.latestValue

    if performTest:
        assertMaximumErrorAngle(Quaternion(),
            filter.rotation.latestValue, 2)

    return estimatedQuaternions

def testOrientationFilters():
    for i in range(3):
        t = np.arange(0, 10, dt)
        trajectory = RandomRotationTrajectory(t)
        for imp in getImplementations(orientation,
                orientation.OrientationFilter):
            params = filterParameters.get(imp, {})
            filter = imp(
                    initialTime=trajectory.startTime,
                    initialRotation=trajectory.rotation(trajectory.startTime),
                    initialRotationalVelocity=trajectory.rotation(trajectory.startTime)
                        .rotateFrame(trajectory.rotationalVelocity(trajectory.startTime)),
                    **params)
            yield checkOrientationFilter, filter, trajectory

def testStaticConvergence():
    for filterClass in getImplementations(orientation,
            orientation.OrientationFilter):
        if filterClass not in [orientation.GyroIntegrator,
                orientation.DistLinAccelCF]:
            yield checkStaticConvergence, filterClass

def testDynamicConvergence():
    for i in range(3):
        t = np.arange(0, 10, dt)
        trajectory = RandomRotationTrajectory(t)
        for imp in getImplementations(orientation,
                orientation.OrientationFilter):
            params = filterParameters.get(imp, {})
            filter = imp(initialRotation=Quaternion(),
                    initialTime=trajectory.startTime,
                    initialRotationalVelocity=trajectory.rotation(trajectory.startTime)
                        .rotateFrame(trajectory.rotationalVelocity(trajectory.startTime)),
                    **params)
            if imp not in [orientation.GyroIntegrator, orientation.DistLinAccelCF]:
                yield checkDynamicConvergence, filter, trajectory
