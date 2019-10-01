"""
Algorithms for reconstructing the posture of a body model.
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

from abc import ABC, abstractmethod
from imusim.trajectories.rigid_body import SampledBodyModel
from imusim.maths.quaternions import Quaternion
from imusim.platforms.radios import RadioPacket


class PostureReconstructor(ABC):
    """
    Base class for posture reconstruction algorithms.

    A posture estimator takes data from IMUs on a jointed rigid body and
    updates the joint rotations of a L{SampledBodyModel}.
    """
    def __init__(self, bodyModel, initialJointRotations=()):
        """
        Initialise posture reconstructor.

        @param bodyModel: The L{SampledBodyModel} to update.
        @param initialJointRotations: List of initial joint rotations
            in depth first pre-order, as L{Quaternion} objects.
        """
        self._bodyModel = bodyModel
        for joint, rotation in zip(bodyModel.joints, initialJointRotations):
            joint.rotationKeyFrames.add(0, rotation)

        self._lastCallTime = None

    def __call__(self, data, t):
        """
        Update the posture of the model given the current model state
        and additional data.

        @param data: Data friom sensors attached to the body model, as a list
            of L{dict} objects. E.g. these may be L{RadioPacket} objects with
            data transmitted by each IMU.
        @param t: The time at which this data was valid.
        """
        dt = t - self._lastCallTime if self._lastCallTime is not None else None
        data = dict((d['jointName'], d) for d in data)
        for joint in self._bodyModel.joints:
            self._update(joint, data.get(joint.name, None), dt, t)

    @abstractmethod
    def _update(self, joint, data, dt, t):
        """
        Internal method to be implemented by subclasses.

        Update the rotation of a joint given data from a sensor attached
        to that joint.

        @param joint:
        @param data: Data from sensors attached to the body model, as a
            L{dict} object. E.g. this may be a L{RadioPacket} with data
            transmitted by an IMU.
        @param dt: The time at which to compute the update (float).
        @param dt: Time elapsed since the last update, or None on first update.
        """
        pass


class SimpleForwardKinematics(PostureReconstructor):
    """
    Posture reconstructor using joint orientations directly.
    """

    def __init__(self, bodyModel):
        PostureReconstructor.__init__(self, bodyModel)

    def _update(self, joint, data, dt, t):
        if data is not None and 'rotation' in data:
            joint.rotationKeyFrames.add(t, data['rotation'])


class InheritedForwardKinematics(PostureReconstructor):
    """
    Posture reconstructor using joint orientations with inheritance fallback.

    Posture reconstruction is performed by applying joint orientations directly
    to the body model, as in SimpleForwardKinematics, however if no data is
    available for a joint it inherits the orientation of its parent joint.
    """

    def __init__(self, bodyModel):
        PostureReconstructor.__init__(self, bodyModel)

    def _update(self, joint, data, dt, t):
        if data is not None and 'rotation' in data:
            joint.rotationKeyFrames.add(t, data['rotation'])
        elif joint.hasParent:
            joint.rotationKeyFrames.add(t, joint.parent.rotation(t))

