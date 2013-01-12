"""
Classes to represent marker-based capture data.
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

from imusim.trajectories.base import AbstractTrajectory
from imusim.trajectories.sampled import SampledPositionTrajectory
from imusim.trajectories.sampled import SampledRotationTrajectory
from imusim.trajectories.splined import SplinedPositionTrajectory
from imusim.trajectories.splined import SplinedRotationTrajectory
import numpy as np

class MarkerCapture(object):
    """
    Marker trajectory data captured synchronously for one or more markers.

    @ivar sampleTimes: Sequence of times at which samples were taken.
    """

    def __init__(self):
        """
        Initialise capture.
        """
        self._markers = dict()

    def _addMarker(self, marker):
        self._markers[marker.id] = marker

    @property
    def markers(self):
        """
        List of the markers in this capture.
        """
        return self._markers.values()

    def marker(self, id):
        """
        Get a marker by identifier.

        @return: The L{Marker} with the given identifier.
        """
        return self._markers[id]

class Marker(AbstractTrajectory):
    """
    A motion capture marker of any type.
    """

    def __init__(self, capture, id):
        """
        Initialise marker.

        @param capture: Parent L{MarkerCapture}.
        @param id: Marker identifier.
        """
        self.capture = capture
        self.id = id
        self.capture._addMarker(self)

class Marker3DOF(SampledPositionTrajectory, Marker):
    """
    A capture marker with 3DOF position samples.
    """
    def __init__(self, capture, id):
        Marker.__init__(self, capture, id)
        SampledPositionTrajectory.__init__(self)

class Marker6DOF(SampledRotationTrajectory, Marker3DOF):
    """
    A capture marker with 6DOF position & orientation samples.
    """
    def __init__(self, capture, id):
        Marker3DOF.__init__(self, capture, id)
        SampledRotationTrajectory.__init__(self)

class SplinedMarkerCapture(MarkerCapture):
    """
    Wrapper to obtain splined versions of all marker trajectories in a capture.

    @param sampled: L{MarkerCapture} with sampled marker data.
    @param kwargs: Keyword arguments to pass to spline trajectory constructors.
    """
    def __init__(self, sampled, **kwargs):
        self.sampled = sampled
        self._markers = dict()
        for marker in sampled.markers:
            if isinstance(marker, Marker6DOF):
                markerClass = SplinedMarker6DOF
            else:
                markerClass = SplinedMarker3DOF
            self._markers[marker.id] = markerClass(self, marker, **kwargs)

class SplinedMarker3DOF(SplinedPositionTrajectory, Marker):
    """
    A capture marker with splined 3DOF position information.
    """
    def __init__(self, capture, sampled, **kwargs):
        """
        Initialise splined marker trajectory.

        @param capture: Parent L{SplinedMarkerCapture}.
        @param sampled: L{Marker3DOF} with sampled data.
        @param kwargs: Keyword arguments for spline trajectory constructors.
        """
        Marker.__init__(self, capture, sampled.id)
        SplinedPositionTrajectory.__init__(self, sampled, **kwargs)

class SplinedMarker6DOF(SplinedRotationTrajectory, SplinedMarker3DOF):
    """
    A capture marker with splined 6DOF position and rotation information.
    """
    def __init__(self, capture, sampled, **kwargs):
        """
        Initialise splined marker trajectory.

        @param capture: Parent L{SplinedMarkerCapture}.
        @param sampled: L{Marker6DOF} with sampled data.
        @param kwargs: Keyword arguments for spline trajectory constructors.
        """
        SplinedMarker3DOF.__init__(self, capture, sampled, **kwargs)
        SplinedRotationTrajectory.__init__(self, sampled, **kwargs)

    @property
    def startTime(self):
        return max(self._positionStartTime, self._rotationStartTime)

    @property
    def endTime(self):
        return min(self._positionEndTime, self._rotationEndTime)
