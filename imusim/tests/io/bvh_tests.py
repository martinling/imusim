"""
Tests for BVH IO
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

from imusim.io import loadBVHFile, saveBVHFile
from imusim.maths.vectors import vector
from imusim.maths.quaternions import Quaternion
from numpy.testing import assert_almost_equal
from imusim.testing.quaternions import assertQuaternionAlmostEqual
from imusim.maths.transforms import convertCGtoNED
from nose.tools import raises
import tempfile

def testBVHInput():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
    JOINT j1
    {
        OFFSET 10.0 0.0 0.0
        CHANNELS 3 Zrotation Xrotation Yrotation
        End Site
        {
            OFFSET 0.0 10.0 0.0
        }
    }
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0 0.0 0.0 0.0 90.0 0.0 0.0
    """
    testFile = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    with testFile:
        testFile.write(data)
        testFile.flush()
        model = loadBVHFile(testFile.name)

        assert len(list(model.points)) == 3

        assert model.name == 'root'
        assert len(model.channels) == 6
        assert len(model.children) == 1
        assert_almost_equal(model.position(0),vector(0,0,0))
        assertQuaternionAlmostEqual(model.rotation(0),
               convertCGtoNED(Quaternion(1,0,0,0)))

        j1 = model.getJoint('j1')
        assert j1.parent is model
        assert len(j1.channels) == 3
        assert len(j1.children) == 1
        assert_almost_equal(j1.positionOffset, vector(10,0,0))
        assert_almost_equal(j1.position(0), convertCGtoNED(vector(10,0,0)))
        assertQuaternionAlmostEqual(j1.rotation(0),
               convertCGtoNED(Quaternion.fromEuler((90, 0, 0), order='zxy')))

        j1end = model.getPoint('j1_end')
        assert j1end.parent is j1
        assert len(j1end.channels) == 0
        assert not j1end.isJoint
        assert_almost_equal(j1end.positionOffset, vector(0,10,0))
        assert_almost_equal(j1end.position(0), vector(0,0,0))

def testBVH_IO():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
    JOINT j1
    {
        OFFSET 10.0 0.0 0.0
        CHANNELS 3 Zrotation Xrotation Yrotation
        End Site
        {
            OFFSET 0.0 10.0 0.0
        }
    }
}
MOTION
Frames: 5
Frame Time: 0.1
0.0 0.0 0.0 0.0 0.0 0.0 90.0 0.0 0.0
0.0 1.0 0.0 0.0 0.0 0.0 80.0 0.0 5.0
0.0 2.0 0.0 0.0 0.0 0.0 70.0 0.0 10.0
0.0 3.0 0.0 0.0 0.0 0.0 60.0 0.0 15.0
0.0 4.0 0.0 0.0 0.0 0.0 50.0 0.0 20.0
    """
    testFile = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    exportFile = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    with testFile:
        testFile.write(data)
        testFile.flush()
        model = loadBVHFile(testFile.name)
        saveBVHFile(model, exportFile.name, 0.1)
        exportedModel = loadBVHFile(exportFile.name)

        assert len(list(model.points)) == len(list(exportedModel.points))
        for orig, exported in zip(model.points, exportedModel.points):
            assert orig.name == exported.name
            assert_almost_equal(orig.positionOffset, exported.positionOffset)
            assert_almost_equal(orig.position(0), exported.position(0))
            if orig.isJoint:
                assertQuaternionAlmostEqual(orig.rotation(0),
                        exported.rotation(0))

def runSyntaxTest(data):
    testFile = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    with testFile:
        testFile.write(data)
        testFile.flush()
        model = loadBVHFile(testFile.name)

@raises(SyntaxError)
def testUnknownChannel():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6 Xpos Yposition Zposition Zrotation Xrotation Yrotation
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0 0.0 0.0 0.0
    """
    runSyntaxTest(data)

@raises(SyntaxError)
def testNoRootPosition():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 3 Zrotation Xrotation Yrotation
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0
    """
    runSyntaxTest(data)

@raises(SyntaxError)
def testNoJointRotation():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
    JOINT j1
    {
        OFFSET 10.0 0.0 0.0
        CHANNELS 0
        End Site
        {
            OFFSET 0.0 10.0 0.0
        }
    }
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0 0.0 0.0 0.0
    """
    runSyntaxTest(data)

@raises(SyntaxError)
def testInsufficientData():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0
    """
    runSyntaxTest(data)

@raises(SyntaxError)
def testInvalidToken():
    data = r"""Broken"""
    runSyntaxTest(data)

@raises(SyntaxError)
def testFloatTokenCheck():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET x 0.0 0.0
    CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0
    """
    runSyntaxTest(data)

@raises(SyntaxError)
def testIntTokenCheck():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6.0 Xposition Yposition Zposition Zrotation Xrotation Yrotation
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0
    """
    runSyntaxTest(data)


@raises(SyntaxError)
def testUnkownToken():
    data = r"""HIERARCHY
ROOT root
{
    OFFSET 0.0 0.0 0.0
    CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
    JOINT j1
    {
        OFFSET 10.0 0.0 0.0
        CHANNELS 3 Zrotation Xrotation Yrotation
        BAD_TOKEN
        End Site
        {
            OFFSET 0.0 10.0 0.0
        }
    }
}
MOTION
Frames: 1
Frame Time: 0.1
0.0 0.0 0.0 0.0 0.0 0.0 90.0 0.0 0.0
    """
    runSyntaxTest(data)
