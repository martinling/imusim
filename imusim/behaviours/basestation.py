"""
Behaviours for basestation devices.
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

from imusim.algorithms.posture_reconstruction import \
    PostureReconstructor, SimpleForwardKinematics
from imusim.algorithms.position import PositionEstimator, ConstantPosition
from imusim.trajectories.rigid_body import SampledBodyModel
import numpy as np

class BodyModelReconstructor(object):
    """
    Reconstructs body model movements from individual joint data.

    This behaviour combines a L{PostureReconstructor} and
    a L{PositionEstimator} to update a L{SampledBodyModel}.

    @ivar bodyModel: L{SampledBodyModel} being updated.
    """

    def __init__(self, bodyModel, initialPosition=np.zeros((3,1)),
            postureReconstructor=SimpleForwardKinematics,
            positionEstimator=ConstantPosition, initialTime=0):
        """
        Initialise behaviour.

        @param bodyModel: L{SampledBodyModel} to update.
        @param postureReconstructor: L{PostureReconstructor} subclass to use.
        @param positionEstimator: L{PositionEstimator} subclass to use.
        @param initialTime: Initial time.
        """
        self.bodyModel = bodyModel
        self._postureEstimator = postureReconstructor(self.bodyModel)
        self._positionEstimator = positionEstimator(self.bodyModel,
                initialPosition = initialPosition)
        self._time = initialTime

    def handleData(self, data, dt):
        """
        Update the body model with new data.

        @param data: A L{dict} containing data, passed directly to both the
            posture and position estimators.
        @param dt: Time elapsed since last call to this method or __init__.
        """
        self._time += dt
        self._postureEstimator(data, self._time)
        self._positionEstimator(data, self._time)
