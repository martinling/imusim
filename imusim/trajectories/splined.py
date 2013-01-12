"""
Continuous splined trajectories from sampled trajectories.
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
from imusim.trajectories.sampled import SampledPositionTrajectory, SampledRotationTrajectory
from imusim.trajectories.sampled import SampledTrajectory
from imusim.maths.vector_splines import PartialInputVectorSpline
from imusim.maths.quat_splines import PartialInputQuaternionBSpline
from imusim.utilities.caching import CacheLastValue
from imusim.utilities.time_series import TimeSeries

class SplinedPositionTrajectory(SampledPositionTrajectory):
    """
    Trajectory with position obtained by splining a sampled trajectory.
    """

    def __init__(self, sampled, positionStdDev=0.001, positionSplineOrder=5, **kwargs):
        """
        Initialise trajectory.

        @param sampled: position samples (L{SampledPositionTrajectory})
        @param positionStdDev: standard deviation of position samples (float)
        @param positionSplineOrder: order of B-spline used for interpolation
        """
        self.sampled = sampled
        self._positionSpline = PartialInputVectorSpline(
            sampled.positionKeyFrames.timestamps,
            sampled.positionKeyFrames.values,
            stddev=positionStdDev, order=positionSplineOrder)

    @property
    def _positionStartTime(self):
        return self._positionSpline.validFrom

    @property
    def _positionEndTime(self):
        return self._positionSpline.validTo

    @property
    def positionKeyFrames(self):
        t = self.sampled.positionKeyFrames.timestamps
        t = t[(t >= self.startTime) & (t <= self.endTime)]
        p = self.position(t)
        return TimeSeries(t, p)

    @CacheLastValue()
    def position(self, t):
        return self._positionSpline(t, 0)

    @CacheLastValue()
    def velocity(self, t):
        return self._positionSpline(t, 1)

    @CacheLastValue()
    def acceleration(self, t):
        return self._positionSpline(t, 2)

    @property
    def startTime(self):
        return self._positionStartTime

    @property
    def endTime(self):
        return self._positionEndTime

class SplinedRotationTrajectory(SampledRotationTrajectory):
    """Trajectory with rotation obtained by splining a sampled trajectory."""

    def __init__(self, sampled, smoothRotations=True, rotationStdDev=0.001, **kwargs):
        """
        Initialise trajectory.

        @param sampled: rotation samples (L{SampledPositionTrajectory})
        @param smoothRotations: `True` to apply smoothing to rotation prior
            to splining (Bool)
        @param rotationStdDev: standard deviation of rotation samples (float)
        """
        self.sampled = sampled
        keyFrames = sampled.rotationKeyFrames.values.unflipped()
        if smoothRotations:
            keyFrames = keyFrames.smoothed(stddev=rotationStdDev)
        self._rotationSpline = PartialInputQuaternionBSpline(
            sampled.rotationKeyFrames.timestamps, keyFrames)

    @property
    def _rotationStartTime(self):
        return self._rotationSpline.validFrom

    @property
    def _rotationEndTime(self):
        return self._rotationSpline.validTo

    @property
    def rotationKeyFrames(self):
        t = self.sampled.rotationKeyFrames.timestamps
        t = t[(t >= self.startTime) & (t <= self.endTime)]
        r = self.rotation(t)
        return TimeSeries(t, r)

    @CacheLastValue()
    def rotation(self, t):
        return self._rotationSpline(t)[0]

    @CacheLastValue()
    def rotationalVelocity(self, t):
        r, omega, alpha = self._rotationSpline(t)
        return r.rotateVector(omega)

    @CacheLastValue()
    def rotationalAcceleration(self, t):
        r, omega, alpha = self._rotationSpline(t)
        return r.rotateVector(alpha)

    @property
    def startTime(self):
        return self._rotationStartTime

    @property
    def endTime(self):
        return self._rotationEndTime

class SplinedTrajectory(SplinedPositionTrajectory, SplinedRotationTrajectory):
    """
    Represents a splined position and rotation trajectory.
    """

    def __init__(self, sampled, **kwargs):
        """
        Initialise trajectory.

        @param sampled: L{SampledTrajectory}

        Keyword arguments are passed to L{SplinedPositionTrajectory} and
        L{SplinedRotationTrajectory} constructors.
        """
        SplinedPositionTrajectory.__init__(self, sampled, **kwargs)
        SplinedRotationTrajectory.__init__(self, sampled, **kwargs)

    @staticmethod
    def fromArrays(time, position, rotation):
        """
        Construct a L{SplinedTrajectory} from time, position and rotation arrays.

        @param time: sequence of sample times.
        @param position: sequence of position samples.
        @param rotation: sequence of rotation samples.
        """
        return SplinedTrajectory(SampledTrajectory.fromArrays(time, position, rotation))

    @property
    def startTime(self):
        return max(self._positionStartTime, self._rotationStartTime)

    @property
    def endTime(self):
        return min(self._positionEndTime, self._rotationEndTime)
