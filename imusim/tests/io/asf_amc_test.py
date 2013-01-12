"""
Tests for ASF/AMC input
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

from imusim.io.asf_amc import loadASFFile
from imusim.io.bvh import loadBVHFile
from os import path
from numpy.testing import assert_almost_equal
from imusim.testing.quaternions import assertMaximumRMSErrorAngle

def testASFInput():
    dir = path.dirname(__file__)
    asfFile = path.join(dir, "16.asf")
    amcFile = path.join(dir, "16_15.amc")
    bvhFile = path.join(dir, "16_15.bvh")

    asfModel = loadASFFile(asfFile, amcFile, scaleFactor=2.54/100,
            framePeriod=1.0/120)
    bvhModel = loadBVHFile(bvhFile, 1.0/100)

    assert set(j.name for j in asfModel) == set(j.name for j in bvhModel)
    for asfNode in asfModel:
        if asfNode == asfModel:
            # skip the root as position offsets are defined in different ways
            continue
        bvhNode = bvhModel.getPoint(asfNode.name)
        assert asfNode.parent.name == bvhNode.parent.name
        assert_almost_equal(asfNode.positionOffset, bvhNode.positionOffset,
                err_msg='Wrong offset for joint %s'%asfNode.name, decimal=6)

    assert_almost_equal(asfModel.positionKeyFrames.values,
            bvhModel.positionKeyFrames.values, decimal=4,
            err_msg="Model position trajectories do not match.")

    for asfJoint in asfModel.joints:
        bvhJoint = bvhModel.getJoint(asfJoint.name)
        assertMaximumRMSErrorAngle(asfJoint.rotationKeyFrames.values,
                bvhJoint.rotationKeyFrames.values, maxRMSE=5,
                errMsg="Joint %s"%asfJoint.name)

