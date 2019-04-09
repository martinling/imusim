"""
Matrix maths tests.
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

from imusim.maths.matrices import matrixToEuler, matrixFromEuler
from imusim.maths.quaternions import QuaternionFromEuler
from imusim.testing.random_data import randomRotation
from imusim.testing.quaternions import assertQuaternionAlmostEqual
import numpy as np
from numpy import testing
from nose.tools import raises

def checkMatrixToEuler(order, degrees):
    q1 = randomRotation()
    m = q1.toMatrix()
    e = matrixToEuler(m, order, degrees)
    q2 = QuaternionFromEuler(e, order, degrees)
    assertQuaternionAlmostEqual(q1,q2)

def testRandomisedMatrixToEuler():
    for degrees in [True, False]:
        for order in ['zyx','zxy']:
            for i in range(10):
                yield checkMatrixToEuler, order, degrees

MATRIX_TO_EULER_TESTS = [
        (np.matrix("1,0,0; 0,1,0; 0,0,1").T, (0,0,0), (0,0,0)),
        (np.matrix("1,0,0; 0,0,1; 0,-1,0").T, (0,0,90), (0,90,0)),
        (np.matrix("1,0,0; 0,0,-1; 0,1,0").T, (0,0,-90), (0,-90,0)),
        (np.matrix("0,0,-1; 0,1,0; 1,0,0").T, (0,90,0), (0,0,90)),
        (np.matrix("0,0,1; 0,1,0; -1,0,0").T, (0,-90,0), (0,0,-90)),
        (np.matrix("0,1,0; -1,0,0; 0,0,1").T, (90,0,0), (90,0,0)),
        (np.matrix("0,-1,0; 1,0,0; 0,0,1").T, (-90,0,0), (-90,0,0))
        ]

def checkZYX(matrix, euler):
    testing.assert_almost_equal(matrixToEuler(matrix, order='zyx',
        inDegrees=True),euler)
    testing.assert_almost_equal(matrixFromEuler(euler, order='zyx'),
            matrix)

def checkZXY(matrix, euler):
    testing.assert_almost_equal(matrixToEuler(matrix, order='zxy',
        inDegrees=True), euler)
    testing.assert_almost_equal(matrixFromEuler(euler, order='zxy'),
            matrix)

def testMatrixToEuler():
    for matrix, zyxAngles, zxyAngles in MATRIX_TO_EULER_TESTS:
        yield checkZXY, matrix, zxyAngles
        yield checkZYX, matrix, zyxAngles

@raises(NotImplementedError)
def testInvalidRotationOrder():
    matrixToEuler(np.eye(3), order='invalid')


