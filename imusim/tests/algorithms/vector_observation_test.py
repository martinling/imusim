"""
Tests for vector observation algorithms.
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
import itertools
import math

from imusim.algorithms import vector_observation
from imusim.maths.quaternions import Quaternion
from imusim.testing.quaternions import assertQuaternionAlmostEqual
from imusim.testing.inspection import getImplementations

angles = [45*p for p in range(8)]

def checkVectorObservation(vectorObservationMethod, eulerAngles, inclination):
    if issubclass(vectorObservationMethod,
            vector_observation.LeastSquaresOptimalVectorObservation):
        vo = vectorObservationMethod(inclinationAngle=inclination)
    else:
        vo = vectorObservationMethod()
    q = Quaternion.fromEuler(eulerAngles)
    inclination = math.radians(inclination)
    m = np.array([[math.cos(inclination)],[0],
        [math.sin(inclination)]])
    g = np.array([[0,0,1]],dtype=float).T
    g = q.rotateFrame(g)
    m = q.rotateFrame(m)
    assertQuaternionAlmostEqual(q,vo(g,m))

def testVectorObservationMethods():
    for method in getImplementations(vector_observation,
            vector_observation.VectorObservation):
        for eulerAngles in itertools.combinations(angles*3,3):
            for inclination in (0,45,80):
                yield checkVectorObservation, method, eulerAngles, \
                        inclination



