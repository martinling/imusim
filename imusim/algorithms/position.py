"""
Algorithms for tracking the position of a body model.
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
import numpy as np
import scipy.stats
import math
from imusim.maths import integrators
from imusim.maths.kalman import KalmanFilter
from imusim.utilities.time_series import TimeSeries
from imusim.utilities.documentation import prepend_method_doc


class PositionEstimator(ABC):
    """
    Base class for position estimation algorithms.

    A position estimator takes data from IMUs on a jointed rigid body and
    updates the root position of a L{SampledBodyModel}.
    """
    def __init__(self, model, initialTime=0, initialPosition=np.zeros((3,1))):
        """
        Initialise position estimator.

        @param model: L{SampledBodyModel} to update.
        @param initialTime: The time at which the the initialPosition is
            sampled
        @param initialPosition: Initial position vector (3x1 L{np.ndarray}).
        """
        self._model = model
        self._model.positionKeyFrames.add(initialTime, initialPosition)

    def __call__(self, data, t):
        """
        Update the position of the model given the current model state
        and additional data.

        @param data: Data from IMUs attached to the body , as a list of
            L{dict} objects. E.g. these may be L{RadioPacket} objects with
            data trasmitted by each IMU.
        @param t: The time at which to compute the update (float).
        """
        dt = t - self._model.positionKeyFrames.latestTime
        position, covariance = self._update(data, dt, t)
        self._model.positionKeyFrames.add(t, position, covariance)

    @abstractmethod
    def _update(self, data, dt, t):
        """
        Internal method to be implemented by subclasses.

        @param data: Data from sensors attached to the body model, as a list
            of L{dict} objects. E.g. these may be L{RadioPacket} objects
            with data trasmitted by each IMU.
        @param dt: Time elapsed since the last update.
        @param t: The time for which to compute the update (float).
        """
        pass

    def getParameter(self, data, jointName, parameter, default=None):
        """
        Get a named parameter for a joint from the passed data tuple.

        @param data: List or tuple of data objects from sensors attached
            to the body model. Data objects are assumed to have a 'jointName'
            attribute with the name of the associated joint and parameters as
            other named attributes e.g. 'acceleration'.
        @param jointName: The joint name to get data for (string).
        @param parameter: The name of the parameter to get data for (string).
        @param default: Default parameter value to return if the parameter is
            not found in the passed data.

        @return: The requested parameter value, or default value if the
            parameter is not found in the passed data.
        """
        for item in data:
            if item.get('jointName') == jointName:
                return item.get(parameter, default)


class ConstantPosition(PositionEstimator):
    """
    Trivial position estimator that gives a constant position.
    """

    def __init__(self, model, initialTime=0, initialPosition=np.zeros((3,1)),
            **kwargs):
        PositionEstimator.__init__(self, model, initialTime, initialPosition)
        self._position = initialPosition.copy()

    def _update(self, data, dt, t):
        return self._position, None

class RootAccelerationIntegrator(PositionEstimator):
    """
    Estimates position based on integration of root joint acceleration.
    """

    @prepend_method_doc(PositionEstimator)
    def __init__(self, model, initialTime=0, initialPosition=np.zeros((3,1)),
            initialVelocity=np.zeros((3,1)),
            integrationMethod=integrators.TrapeziumRule, **kwargs):
        """
        @param initialVelocity: Initial velocity estimate (3x1 L{np.ndarray})
        @param integrationMethod: L{Integrator} to use for position integration.
        """
        PositionEstimator.__init__(self, model, initialTime, initialPosition)
        self._integrator = integrators.DoubleIntegrator(initialPosition,
                initialVelocity, integrationMethod)

    def _update(self, data, dt, t):
        acceleration = self.getParameter(data, self._model.name,
                'linearAcceleration', default=np.zeros((3,1)))

        return self._integrator(acceleration, dt), None

class LowestPoint(PositionEstimator):
    """
    Tracking based on assumption that lowest point is in static ground contact.
    """

    def __init__(self, model, initialTime=0, initialPosition=np.zeros((3,1)),
            **kwargs):
        PositionEstimator.__init__(self, model, initialTime, initialPosition)
        self._anchor = self._lowestPoint(0)
        self._anchorPosition = self._anchor.position(0)

    def _lowestPoint(self, t):
        key = lambda point: point.position(t)[2]
        return sorted([p for p in self._model.points if
            not np.any(np.isnan(p.position(t)))], key=key)[-1]

    def _update(self, data, dt, t):
        lowest = self._lowestPoint(t)
        if self._anchor is not lowest:
            self._anchorPosition += (lowest.position(t) -
                self._anchor.position(t))
            self._anchorPosition[2] = 0
            self._anchor = lowest

        return self._anchorPosition + (self._model.position(t) -
                self._anchor.position(t)), None

class ContactTrackingKalmanFilter(PositionEstimator):
    """
    Estimates position based on root acceleration and estimated ground contact.
    """

    @prepend_method_doc(PositionEstimator)
    def __init__(self, model, initialTime=0,
            initialPosition=np.zeros((3,1)),
            initialVelocity=np.zeros((3,1)),
            accelerationVar=0.2, velocityVar=0.5,
            rejectionProbabilityThreshold=0.4,
            **kwargs):
        """
        @param initialVelocity: Initial velocity estimate (3x1 L{np.ndarray})
        @param accelerationVar: Acceleration variance (float).
        @param velocityVar: Velocity variance (float).
        @param rejectionProbabilityThreshold: Probability below which a
            velocity update will be rejected (float, 0 <= p <= 1)
        """
        PositionEstimator.__init__(self, model, initialTime, initialPosition)

        self._accelerationSigma = math.sqrt(accelerationVar)
        self._measurementVariance = np.diag([velocityVar]*3)
        self._rejectionThreshold = rejectionProbabilityThreshold

        self._points = list(model.points)
        for point in self._points:
            point.contactProbabilities = TimeSeries()

        self.kalman = KalmanFilter(
                stateTransitionMatrix = self.transitionMatrix(1),
                controlMatrix = self.controlMatrix(1),
                measurementMatrix = np.hstack((np.zeros((3,3)),np.eye(3))),
                state = np.vstack((initialPosition,initialVelocity)),
                stateCovariance = self.processCovariance(1),
                processCovariance = self.processCovariance(1),
                measurementCovariance = self._measurementVariance)

    def _update(self, data, dt, t):
        acceleration = self.getParameter(data, self._model.name,
                'linearAcceleration', default=np.zeros((3,1)))

        self.kalman.stateTransitionMatrix = self.transitionMatrix(dt)
        self.kalman.controlMatrix = self.controlMatrix(dt)
        self.kalman.processCovariance = self.processCovariance(dt)

        self.kalman.predict(acceleration)

        velocity, probability = self.estimateVelocityProbabilities(dt, t)
        if probability > self._rejectionThreshold:
            self.kalman.correct(velocity)

        return self.kalman.state[:3], self.kalman.stateCovariance[:3, :3]

    def estimateVelocityProbabilities(self, dt, t):
        """
        Estimate the possible velocities (and probabilities) of the model.

        Estimates are based on the assumption that a point is stationary
        in world co-ordinates in the interval `[t-dt,t]`

        @return: The velocity with the highest probabilty of not being wrong
        """
        vPredicted = self.kalman.state[3:]
        inv_cov = self.kalman.stateCovariance[3:6,3:6].I

        highestProbability = 0

        for point in self._points:
            vHat = -(point.position(t) - point.position(t-dt)) / dt
            vError = vHat - vPredicted
            mahalanobis_distance = np.sqrt(vError.T * inv_cov * vError)
            probability = 2 * scipy.stats.norm.cdf(-mahalanobis_distance)
            if probability > highestProbability:
                highestProbability = probability
                estimatedVelocity = vHat
            point.contactProbabilities.add(t, probability[0,0])

        return estimatedVelocity, highestProbability

    def transitionMatrix(self, dt):
        A = np.matrix([ [1, 0, 0, dt, 0, 0],
                        [0, 1, 0, 0, dt, 0],
                        [0, 0, 1, 0, 0, dt],
                        [0, 0, 0, 1, 0, 0 ],
                        [0, 0, 0, 0, 1, 0 ],
                        [0, 0, 0, 0, 0, 1 ]])
        return A

    def controlMatrix(self, dt):
        a = 0.5 * dt**2
        b = dt
        B = np.matrix([ [a, 0, 0],
                        [0, a, 0],
                        [0, 0, a],
                        [b, 0, 0],
                        [0, b, 0],
                        [0, 0, b]])
        return B

    def processCovariance(self, dt):
        a = (0.5 * dt**2)**2 * self._accelerationSigma**2
        b = ((0.5 * dt**2) * self._accelerationSigma) * (dt *
                self._accelerationSigma)
        c = dt**2 * self._accelerationSigma**2

        Q = np.matrix([ [a, 0, 0, b, 0, 0],
                        [0, a, 0, 0, b, 0],
                        [0, 0, a, 0, 0, b],
                        [b, 0, 0, c, 0, 0],
                        [0, b, 0, 0, c, 0],
                        [0, 0, b, 0, 0, c]])
        return Q
