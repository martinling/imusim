"""
Tests for position estimation algorithms.
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

from os import path
from imusim.io.bvh import loadBVHFile
from imusim.trajectories.rigid_body import SampledBodyModel, SplinedBodyModel
from imusim.algorithms import position
from imusim.testing.inspection import getImplementations
import numpy as np

referenceModel = SplinedBodyModel(
    loadBVHFile(path.join(path.dirname(__file__), 'test.bvh'), 0.01))
dt = 0.01
sampleTimes = np.arange(referenceModel.startTime, referenceModel.endTime, dt)

def checkAlgorithm(algorithm):
    sampledModel = SampledBodyModel.structureCopy(referenceModel)

    for t in sampleTimes:
        for referenceJoint, sampledJoint in zip(referenceModel.joints,
                sampledModel.joints):
            sampledJoint.rotationKeyFrames.add(t, referenceJoint.rotation(t))

    positionEstimator = algorithm(sampledModel,
            initialTime = referenceModel.startTime,
            initialPosition = referenceModel.position(referenceModel.startTime),
            initialVelocity = referenceModel.velocity(referenceModel.startTime))

    for t in sampleTimes[1:]:
        data = [{'jointName' : referenceModel.name,
                'linearAcceleration' : referenceModel.acceleration(t)}]
        positionEstimator(data, t)

    result = sampledModel.positionKeyFrames.values
    assert np.shape(result) == referenceModel.position(sampleTimes).shape
    assert np.all(np.isfinite(result))

def testPositionAlgorithms():
    for algorithm in getImplementations(position, position.PositionEstimator):
        yield checkAlgorithm, algorithm
