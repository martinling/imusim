"""
Modelling of rigid body models for motion trajectory generation.

Rigid body models are composed of a tree of L{Joint} objects. Each Joint has
an offset from its parent, specified in the parents local co-ordinate frame.
The position of any part of the body model can be calculated using forward
kinematics, given the position of the root joint of the body model and the
rotation of each of the joints.

The rotation of joints is given by a L{RotationTrajectory} specifying the
rotation of the joint local co-ordinate frame relative to the global
co-ordinate frame at any point in time. The position of the root joint of
the body model is specified by a L{PositionTrajectory}.
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

from imusim.trajectories.offset import OffsetTrajectory
from imusim.trajectories.sampled import SampledPositionTrajectory
from imusim.trajectories.sampled import SampledRotationTrajectory
from imusim.trajectories.splined import SplinedPositionTrajectory
from imusim.trajectories.splined import SplinedRotationTrajectory
from imusim.utilities.trees import TreeNode
import numpy as np

class Point(TreeNode):
    """
    A point in a rigid body model.

    @ivar name: Name of this point.
    @ivar parent: The parent L{Joint} in the body model hierarchy.
    @ivar positionOffset: Position offset of this point from its parent in the
        parent's local co-ordinate frame, or the global frame if no parent.
    """

    def __init__(self, parent, name=None, offset=None):
        """
        Construct a point.

        @param parent: Parent L{Joint}, or None for a point in the global
            co-ordinate frame.
        @param name: Name of this point.
        @param offset: Offset of this point in the co-ordinate frame of its
            parent, as a 3x1 L{np.ndarray}. Default is zero offset.
        """
        TreeNode.__init__(self, parent)
        assert isinstance(self, Joint) or isinstance(parent, Joint)
        self.positionOffset = np.zeros((3,1)) if offset is None else offset
        self.name = name

    @property
    def isJoint(self):
        """
        Whether this point is a joint that can rotate.
        """
        return isinstance(self, Joint)

class Joint(Point):
    """
    A joint in a rigid body model, with its own rotating co-ordinate frame.

    @ivar children: List of child points that have this joint as parent.
    """

    def __init__(self, parent, name=None, offset=None):
        """
        Construct a joint.

        @param parent: Parent L{Joint}, or None for the root joint.
        @param name: Name of this joint.
        @param offset: Offset of this joint in the co-ordinate frame of its
            parent, as a 3x1 L{np.ndarray}. Default is zero offset.
        """
        Point.__init__(self, parent, name, offset)
        self.children = []

    @property
    def points(self):
        """
        Iterator for the points of the body model starting at this point.

        @see: L{preorderTraversal}
        """
        return self.preorderTraversal()

    @property
    def pointNames(self):
        """
        List of point names in the tree rooted at this joint.
        """
        return [p.name for p in self.points]

    def getPoint(self, pointName):
        """
        Get a point by name from the point tree rooted at this point.

        @param pointName: The name of the point to return (string).
        @return: The corresponding L{Point} object.
        @raises KeyError: If the named point is not present in the subtree.
        """
        for p in self.points:
            if p.name == pointName:
                return p
        raise KeyError

    @property
    def childJoints(self):
        """
        List of child joints that have this joint as parent.
        """
        return filter(lambda p: p.isJoint, self.children)

    @property
    def hasChildJoints(self):
        """
        Whether this joint has any child joints.
        """
        return len(self.childJoints) != 0

    @property
    def joints(self):
        """
        Iterator for the joints of the body model starting at this joint.
        """
        return self.preorderTraversal(condition=lambda p: p.isJoint)

    @property
    def jointNames(self):
        """
        List of joint names in the body model tree rooted at this joint.
        """
        return [j.name for j in self.joints]

    def getJoint(self, jointName):
        """
        Get a joint by name from the subtree rooted at this joint.

        @param jointName: The name of the joint to return (string).
        @return: The corresponding L{Joint} object.
        @raises KeyError: If the named joint is not present in the subtree.
        """
        for p in self.joints:
            if p.name == jointName:
                return p
        raise KeyError

    @classmethod
    def structureCopy(cls, ref, parent=None):
        """
        Create a deep copy of the body model tree starting from the given joint.
        """
        joint = cls(parent, ref.name, ref.positionOffset)
        for child in ref.children:
            if child.isJoint:
                cls.structureCopy(child, joint)
            else:
                Point(joint, child.name, child.positionOffset)
        return joint

class PointTrajectory(OffsetTrajectory, Point):
    """
    Trajectory followed by a point in a rigid body model.
    """
    def __init__(self, parent, name=None, offset=None):
        """
        Construct point trajectory.

        @param parent: Parent L{Joint} in the body model hierarchy.
        @param name: Name of the point.
        @param offset: Offset of the point in the parent joint
            co-ordinate frame.
        """
        Point.__init__(self, parent, name, offset)
        OffsetTrajectory.__init__(self, parent, offset)

class SampledJoint(SampledRotationTrajectory, OffsetTrajectory, Joint):
    """
    A joint in a rigid body model with sampled rotations.
    """
    def __init__(self, parent, name=None, offset=None):
        """
        Initialise sampled joint.

        @param parent: Parent L{Joint} in the body model hierarchy.
        @param name: Name of the joint.
        @param offset: Offset of the joint in the parent joint
            co-ordinate frame.
        """
        Joint.__init__(self, parent, name, offset)
        OffsetTrajectory.__init__(self, parent, offset)
        SampledRotationTrajectory.__init__(self)

    @classmethod
    def structureCopy(cls, ref, parent=None):
        joint = cls(parent, ref.name, ref.positionOffset)
        for child in ref.children:
            if child.isJoint:
                cls.structureCopy(child, joint)
            else:
                PointTrajectory(joint, child.name, child.positionOffset)
        return joint

class SampledBodyModel(SampledPositionTrajectory, SampledJoint):
    """A body model with sampled positions and rotations."""

    def __init__(self, name=None):
        """
        Initialise body model.

        @param name: The name of the root L{Joint} of the body model
        """
        SampledJoint.__init__(self, None, name)
        SampledPositionTrajectory.__init__(self)

    @classmethod
    def structureCopy(cls, ref, parent=None):
        model = cls(ref.name)
        for child in ref.children:
            if child.isJoint:
                SampledJoint.structureCopy(child, model)
            else:
                PointTrajectory(model, child.name, child.positionOffset)
        return model

class SplinedJoint(SplinedRotationTrajectory, OffsetTrajectory, Joint):
    """A joint with rotations splined from sampled data."""

    def __init__(self, parent, sampled, **kwargs):
        """
        Initialise splined joint.

        @param sampled: the L{SampledJoint} from which to generate splined
            trajectories
        @param parent: the parent joint in the body model hierarchy
        """
        Joint.__init__(self, parent, sampled.name, sampled.positionOffset)
        OffsetTrajectory.__init__(self, parent, sampled.positionOffset)
        SplinedRotationTrajectory.__init__(self, sampled, **kwargs)
        self.children = [
                (SplinedJoint(self, child, **kwargs) if child.isJoint
                else PointTrajectory(self, child.name, child.positionOffset))
                for child in sampled.children]

    @classmethod
    def structureCopy(cls, ref):
        raise NotImplementedError

    @property
    def _rotationStartTime(self):
        return max([SplinedRotationTrajectory._rotationStartTime.fget(self)] +
            [c._rotationStartTime for c in self.childJoints])

    @property
    def _rotationEndTime(self):
        return min([SplinedRotationTrajectory._rotationEndTime.fget(self)] +
            [c._rotationEndTime for c in self.childJoints])

    @property
    def startTime(self):
        return self.root.startTime

    @property
    def endTime(self):
        return self.root.endTime

class SplinedBodyModel(SplinedPositionTrajectory, SplinedJoint):
    """A body model with positions and rotations splined from sampled data."""

    def __init__(self, sampled, **kwargs):
        """
        Initialise splined body model.

        @param sampled: The L{SampledBodyModel} from which to generated
            splined trajectories
        """
        SplinedJoint.__init__(self, None, sampled, **kwargs)
        SplinedPositionTrajectory.__init__(self, sampled, **kwargs)

    @property
    def startTime(self):
        return max(self._positionStartTime, self._rotationStartTime)

    @property
    def endTime(self):
        return min(self._positionEndTime, self._rotationEndTime)
