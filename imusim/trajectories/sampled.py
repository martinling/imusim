"""
Classes representing sampled trajectories.

Sampled trajectories do not support evaluation of derivative values such as
velocities or accelerations. To evaulate these properties a splined trajectory
should be constructed. See L{splined}
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

from .base import PositionTrajectory, RotationTrajectory
from imusim.maths.quaternions import Quaternion, QuaternionArray
from imusim.utilities.caching import CacheLastValue
from imusim.utilities.time_series import TimeSeries
import imusim.maths.vectors as vectors
import numpy as np

class SampledPositionTrajectory(PositionTrajectory):
    """
    Represents a sampled position trajectory.

    @ivar positionKeyFrames: A L{TimeSeries} of position key frames.
    """

    def __init__(self, keyFrames=None):
        """
        Initialise trajectory.

        @param keyFrames: A L{TimeSeries} of position key frame samples.
        """
        self.positionKeyFrames = TimeSeries() if keyFrames is None else keyFrames

    @CacheLastValue()
    def position(self, t):
        if len(self.positionKeyFrames) == 0:
            if np.ndim(t) == 0:
                return vectors.nan()
            else:
                return vectors.nan(len(t))
        else:
            return self.positionKeyFrames(t)

    def velocity(self, t):
        """
        Not implemented for sampled trajectories.

        See L{imusim.trajectories.splined} for support.
        """
        raise NotImplementedError("Derivative not available from sampled trajectory. " \
            + "Create a splined trajectory to obtain derivatives.")

    def acceleration(self, t):
        """
        Not implemented for sampled trajectories.

        See L{imusim.trajectories.splined} for support.
        """
        raise NotImplementedError("Derivative not available from sampled trajectory. " \
            + "Create a splined trajectory to obtain derivatives.")

    @property
    def startTime(self):
        return self.positionKeyFrames.earliestTime

    @property
    def endTime(self):
        return self.positionKeyFrames.latestTime

class SampledRotationTrajectory(RotationTrajectory):
    """
    Represents a sampled rotation trajectory.

    @ivar rotationKeyFrames: A L{TimeSeries} of rotiaton key frames.
    """

    def __init__(self, keyFrames=None):
        """
        Initialise trajectory.

        @param keyFrames: A L{TimeSeries} of rotation key frame samples.
        """
        self.rotationKeyFrames = TimeSeries() if keyFrames is None else keyFrames

    @CacheLastValue()
    def rotation(self, t):
        if len(self.rotationKeyFrames) == 0:
            if np.ndim(t) == 0:
                return Quaternion.nan()
            else:
                return QuaternionArray.nan(len(t))
        else:
            return self.rotationKeyFrames(t)

    def rotationalVelocity(self, t):
        """
        Not implemented for sampled trajectories.

        See L{splined} for support.
        """
        raise NotImplementedError("Derivative not available from sampled trajectory. " \
            + "Create a splined trajectory to obtain derivatives.")

    def rotationalAcceleration(self, t):
        """
        Not implemented for sampled trajectories.

        See L{splined} for support.
        """
        raise NotImplementedError("Derivative not available from sampled trajectory. " \
            + "Create a splined trajectory to obtain derivatives.")

    @property
    def startTime(self):
        return self.rotationKeyFrames.earliestTime

    @property
    def endTime(self):
        return self.rotationKeyFrames.latestTime

class SampledTrajectory(SampledPositionTrajectory, SampledRotationTrajectory):
    """
    Represents a sampled position and rotation trajectory.
    """

    def __init__(self, positionKeyFrames=None, rotationKeyFrames=None):
        """
        Initialise trajectory.

        @param positionKeyFrames: L{TimeSeries} of position key frames.
        @param rotationKeyFrames: L{TimeSeries} of rotation key frames.
        """
        SampledPositionTrajectory.__init__(self, positionKeyFrames)
        SampledRotationTrajectory.__init__(self, rotationKeyFrames)

    @staticmethod
    def fromArrays(time, position, rotation):
        """
        Construct a L{SampledTrajectory} from time, position and rotation arrays.

        @param time: sequence of sample times.
        @param position: sequence of position samples.
        @param rotation: sequence of rotation samples.
        """
        return SampledTrajectory(
                TimeSeries(time, position),
                TimeSeries(time, rotation))

    @property
    def startTime(self):
        return max(self.positionKeyFrames.earliestTime, self.rotationKeyFrames.earliestTime)

    @property
    def endTime(self):
        return min(self.positionKeyFrames.latestTime, self.rotationKeyFrames.latestTime)
