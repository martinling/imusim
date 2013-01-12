"""
Tests for quaternion spline implementations.
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
from imusim.maths.quat_splines import QuaternionBSpline
from imusim.maths.quat_splines import PartialInputQuaternionBSpline
from imusim.maths.quaternions import Quaternion, QuaternionArray
from imusim.testing.quaternions import assertQuaternionAlmostEqual
from imusim.testing.random_data import randomValidity, invalidate

def geodesicQuaternionPath():
    """
    Generate an array of quaternions following a geodesic path

    Returns
    -------
    t : :class:`~numpy.ndarray`
        timestamps of quaternions
    quaternionPath : :class:`~imusim.maths.quaternions.QuaternionArray`
        array of quaternions following a geodesic path
    """

    t = np.arange(0,10,0.1)
    q = Quaternion.fromEuler((45,0,0))

    return t, QuaternionArray([q**T for T in t])

def checkQuaternionSpline(splineClass, t, qa):
    spline = splineClass(t,qa)
    qs,_,_ = spline(t)
    assertQuaternionAlmostEqual(qa[2:-2], qs[2:-2])

def testQuaternionBSpline():
    t,qa = geodesicQuaternionPath()
    yield checkQuaternionSpline, QuaternionBSpline, t, qa
    for i in range(10):
        validity = randomValidity(t)
        qb = invalidate(qa, validity)
        yield checkQuaternionSpline, PartialInputQuaternionBSpline, t, qb
