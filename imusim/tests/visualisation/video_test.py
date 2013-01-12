"""
Test video generation.
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

from imusim.io.bvh import loadBVHFile
from imusim.trajectories.rigid_body import SampledBodyModel, SplinedBodyModel
from imusim.environment.magnetic_fields import SolenoidMagneticField
from imusim.maths.transforms import AffineTransform
from imusim.algorithms.position import ContactTrackingKalmanFilter
from imusim.visualisation.video import createVideo
from imusim.visualisation.rendering import BodyModelRenderer, \
        VelocityVectorRenderer, ContactProbabilityRenderer, \
        renderSolenoidCoil, autoPositionCamera
import numpy as np
import tempfile
import os
import shutil

def testVideoGeneration():
    samp = loadBVHFile(os.path.join(os.path.dirname(__file__),
        '../system/walk.bvh'), 0.01)
    spl = SplinedBodyModel(samp)

    positionEstimator = ContactTrackingKalmanFilter(
            SampledBodyModel.structureCopy(samp),
            initialTime = spl.startTime,
            initialPosition = spl.position(spl.startTime),
            initialVelocity = spl.velocity(spl.startTime))

    dt = 0.01
    sampleTimes = np.arange(spl.startTime + dt, spl.endTime, dt)

    for t in sampleTimes:
        data = [{'jointName' : p.name,
                'linearAcceleration' : p.acceleration(t)}
            for p in spl.joints]
        positionEstimator(data, t)

    outdir = tempfile.mkdtemp()
    filename = os.path.join(outdir, 'movie.avi')

    autoPositionCamera()
    renderSolenoidCoil(
        SolenoidMagneticField(200,20,0.05,0.2, AffineTransform()))
    createVideo(filename, [BodyModelRenderer(spl),
        VelocityVectorRenderer(spl.getPoint('ltoes_end')),
        ContactProbabilityRenderer(spl)],
        spl.startTime, spl.startTime + 1, 3, (320, 200))
    assert os.path.exists(filename)
    shutil.rmtree(outdir)




