"""
6DOF trajectories from multiple 3DOF marker trajectories.
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

from imusim.trajectories.sampled import SampledPositionTrajectory, SampledRotationTrajectory
from imusim.algorithms.vector_observation import DavenportQ
from imusim.maths.quaternions import QuaternionArray
from imusim.utilities.time_series import TimeSeries
import numpy as np

class MultiMarkerTrajectory(SampledPositionTrajectory, SampledRotationTrajectory):

    def __init__(self, posMarker, refMarkers, refVectors=None, refTime=0):
        self.posMarker = posMarker
        self.refMarkers = refMarkers
        if refVectors is None:
            self.refVectors = [
                refMarker.position(refTime) - posMarker.position(refTime)
                    for refMarker in self.refMarkers]
        else:
            self.refVectors = refVectors
        assert(~np.any(np.isnan(self.refVectors)))

        self.positionKeyFrames = self.posMarker.positionKeyFrames
        t = self.positionKeyFrames.timestamps
        vectorObservation = DavenportQ(self.refVectors)
        measVectors = np.array([ref.position(t) - self.posMarker.position(t)
            for ref in self.refMarkers])
        rotations = QuaternionArray([
            vectorObservation(*vectors).conjugate for vectors in
                np.split(measVectors, len(t), axis=2)])
        self.rotationKeyFrames = TimeSeries(t, rotations)
