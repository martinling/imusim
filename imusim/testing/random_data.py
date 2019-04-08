"""
Functions and classes for generating random data to use in testing.
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
from imusim.maths.quaternions import QuaternionArray, QuaternionFromEuler
from imusim.trajectories.sampled import SampledTrajectory
from imusim.trajectories.splined import SplinedRotationTrajectory
from imusim.trajectories.splined import SplinedTrajectory
from imusim.utilities.time_series import TimeSeries

rng = np.random.RandomState(42)

def randomPosition(range=(-2, 2)):
    """
    A random 3x1 position vector in the given range.
    """
    return rng.uniform(*(list(range) + [(3,1)]))

def randomVelocity(range=(-2, 2)):
    """
    A random 3x1 velocity vector in the given range.
    """
    return rng.uniform(*(list(range) + [(3,1)]))

def randomAngles():
    """
    A random Euler angle sequence in radians.
    """
    return rng.uniform(0, np.pi, (3,1))

def randomRotation():
    """
    A random spatial rotation quaternion.
    """
    return QuaternionFromEuler(randomAngles()[:,0].T, inDegrees=False)

def randomTimeSequence():
    """
    A random sequence of equally spaced time values.
    """
    t = rng.uniform(0, 10)
    l = rng.uniform(0.6, 2)
    dt = rng.uniform(0.01, 0.02)
    return np.arange(t, t+l, dt)

def randomPositionSequence(t):
    """
    A random curving path of position values at given times.
    """
    p = randomPosition()
    v = randomVelocity()
    theta = randomAngles()
    omega = rng.uniform(-10, 10, (3,1))
    a = rng.uniform(-10, 10, (3,1))
    return p + v*t + a*np.sin(omega*t + theta)

def randomRotationSequence(t):
    """
    A random curving path of rotation values at given times.
    """
    theta = randomAngles()
    omega = rng.uniform(-2, 2, (3,1))
    phi = randomAngles()
    omicron = rng.uniform(-10, 10, (3,1))
    a = rng.uniform(-1, 1, (3,1))
    angles = theta + omega*t + a*np.sin(omicron*t + phi)
    return QuaternionArray([
        QuaternionFromEuler(angles[:,i].T, inDegrees=False)
            for i in range(len(t))])

def randomValidity(t):
    """
    A random pattern of data validity at given times.
    """
    validity = np.empty(np.shape(t),dtype=bool)
    valid = (rng.uniform(0,1) > 0.5)
    i = 0
    while i < len(t):
        length = int(rng.uniform(1, 10))
        validity[i:min(i+length,len(validity))] = valid
        valid = ~valid
        i += length
    return validity

def invalidate(data, validity):
    """
    Invalidate data in a given pattern by replacing with NaN values.
    """
    result = data.copy()
    np.atleast_2d(result)[:,~validity] = np.nan
    return result

class RandomRotationTrajectory(SplinedRotationTrajectory):
    """
    A trajectory with a random curving path of rotation values.
    """
    def __init__(self, t):
        sampled = SampledTrajectory(None, TimeSeries(t, randomRotationSequence(t)))
        SplinedRotationTrajectory.__init__(self, sampled)

    def position(self, t):
        return np.zeros((3, len(np.atleast_1d(t))))

    def velocity(self,t):
        return np.zeros((3, len(np.atleast_1d(t))))

    def acceleration(self,t):
        return np.zeros((3, len(np.atleast_1d(t))))

class RandomTrajectory(SplinedTrajectory):
    """
    A trajectory with random curving paths of position and rotation.
    """
    def __init__(self):
        t = randomTimeSequence()
        sampled = SampledTrajectory(
            TimeSeries(t, randomPositionSequence(t)),
            TimeSeries(t, randomRotationSequence(t)))
        SplinedTrajectory.__init__(self, sampled)
