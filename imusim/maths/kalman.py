"""
Standard and unscented Kalman filters.
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

from __future__ import division
import numpy as np
from imusim.maths.transforms import UnscentedTransform
from numpy.linalg import inv

class KalmanFilter(object):
    """
    Implementation of standard Kalman filter

    @ivar stateTransitionMatrix: State transition matrix.
    @ivar controlMatrix: Control input matrix.
    @ivar measurementMatrix: Measurement observation matrix.
    @ivar state: Current state estimate.
    @ivar stateCovariance: State covariance matrix.
    @ivar processCovariance: Process noise covariance matrix.
    @ivar measurementCovariance: Measurement noise covariance matrix.
    """

    def __init__(self, stateTransitionMatrix, controlMatrix, measurementMatrix,
            state, stateCovariance, processCovariance, measurementCovariance):
        """
        Construct Kalman filter.

        @param stateTransitionMatrix: State transition matrix.
        @param controlMatrix: Control input matrix.
        @param measurementMatrix: Measurement observation matrix.
        @param state: Initial state estimate.
        @param stateCovariance: Initial state covariance matrix.
        @param processCovariance: Process noise covariance matrix.
        @param measurementCovariance: Measurement noise covariance matrix.
        """
        self.state = state
        self.stateTransitionMatrix = stateTransitionMatrix
        self.controlMatrix = controlMatrix
        self.measurementMatrix = measurementMatrix
        self.stateCovariance = stateCovariance
        self.processCovariance = processCovariance
        self.measurementCovariance = measurementCovariance

    @property
    def state(self):
        return self._x
    @state.setter
    def state(self, state):
        assert 0 < np.ndim(state) <= 2, "state must be a vector"
        x = state.reshape(-1,1)
        self._states = x.shape[0]
        self._x = np.asmatrix(x)
        self._I = np.eye(self._states)

    @property
    def stateTransitionMatrix(self):
        return self._A
    @stateTransitionMatrix.setter
    def stateTransitionMatrix(self, stateTransitionMatrix):
        assert stateTransitionMatrix.shape == (self._states,self._states), \
            "stateTransitionMatrix must be a %d x %d array/matrix" \
            %(self._states,self._states)
        self._A = np.asmatrix(stateTransitionMatrix)

    @property
    def controlMatrix(self):
        return self._B
    @controlMatrix.setter
    def controlMatrix(self, controlMatrix):
        assert controlMatrix is None or (controlMatrix.ndim == 2 and \
                controlMatrix.shape[0] == self._states), \
                "controlMatrix must be a %d x u.ndim array/matrix, "\
                "where u is the control input vector" % self._states
        self._B = None if controlMatrix is None else np.asmatrix(controlMatrix)

    @property
    def measurementMatrix(self):
        return self._C
    @measurementMatrix.setter
    def measurementMatrix(self, measurementMatrix):
        assert measurementMatrix.ndim == 2 and measurementMatrix.shape[1] == self._states, \
                "measurementMatrix must be a "\
                "z.ndim x %d array/matrix, where z is the observation "\
                "vector" %self._states
        self._C = np.asmatrix(measurementMatrix)
        self._measurements = self._C.shape[0]

    @property
    def stateCovariance(self):
        return self._P
    @stateCovariance.setter
    def stateCovariance(self, stateCovariance):
        assert stateCovariance.shape == (self._states,self._states), "stateCovariance"\
                " must be a %d x %d array/matrix" %(self._states,self._states)
        self._P = np.asmatrix(stateCovariance)

    @property
    def processCovariance(self):
        return self._Q
    @processCovariance.setter
    def processCovariance(self, processCovariance):
        assert processCovariance.shape == (self._states,self._states), \
                "processCovariance must be a %d x %d array/matrix" %(
                self._states,self._states)
        self._Q = np.asmatrix(processCovariance)

    @property
    def measurementCovariance(self):
        return self._R
    @measurementCovariance.setter
    def measurementCovariance(self, measurementCovariance):
        if hasattr(self, '_measurements'):
            assert measurementCovariance.shape == (self._measurements,self._measurements), \
                "measurementCovariance must be a %d x %d array/matrix" %(
                self._measurements,self._measurements)
        self._R = np.asmatrix(measurementCovariance)

    @property
    def gain(self):
        return self._K

    def predict(self, control=None):
        """
        Predict the system state at the next timestep, updating the estimate.

        @param control: Mx1 L{np.ndarray} of control inputs, or None if there
            are no control inputs to the system.
        """
        A,B,P,Q,x = self._A, self._B, self._P, self._Q, self._x

        assert (control is None) == (B is None), \
            "Control input must be given if and only if control matrix is set."

        if B is None:
            B = 0
            u = 0
        else:
            u = np.asanyarray(control).reshape(-1,1)

        self._x = A*x + B*u
        self._P = A*P*A.T + Q

    def innovation(self, measurement):
        """
        Calculate the filter innovation and covariance for a given measurement.

        @param measurement: Nx1 L{np.ndarray} of measurement inputs.

        @return: Nx1 innovation vector, NxN covariance matrix of the
            innovation, and MxN cross covariance from states to measurements.
        """
        C,P,R,x = self._C, self._P, self._R, self._x
        y = np.asanyarray(measurement).reshape(-1,1)

        innovation = y - C*x
        crossCovariance = P*C.T
        innovationCovariance = C*crossCovariance + R

        return innovation, innovationCovariance, crossCovariance

    def correct(self, measurement):
        """
        Correct the filter state given a new measurement.

        @param measurement: Nx1 L{np.ndarray} of measurement inputs.
        """
        y = np.asanyarray(measurement).reshape(-1,1)
        innovation, innovationCovariance, crossCovariance = \
                self.innovation(y)

        self._K = crossCovariance * inv(innovationCovariance)
        self._x = self._x + self._K*innovation
        self._P = (self._I - self._K*self._C) * self._P

    def update(self, measurement, control=None):
        """
        Run filter update given measurement and control inputs.

        @param measurement: Nx1 L{np.ndarray} of measurement inputs.
        @param control: Mx1 L{np.ndarray} of control inputs, or None if there
            are no control inputs to the system.

        """
        self.predict(control)
        self.correct(measurement)

    def normalisedInnovation(self, measurement):
        """
        Obtain the normalised innovation for a given measurement.

        @param measurement: Nx1 L{np.ndarray} of measurement inputs.
        """
        innovation, innovationCovariance, crossCovariance = \
                self.innovation(measurement)
        normalisedInnovation = np.empty_like(innovation)
        for j,i in enumerate(innovation):
            normalisedInnovation[j] = i / np.sqrt(innovationCovariance[j,j])

        return normalisedInnovation

class UnscentedKalmanFilter(KalmanFilter):
    """
    Implementation of the unscented Kalman filter algorithm for additive noise.

    Implementation equations are taken from Wan and Van Der Merwe, "The
    Unscented Kalman Filter", table 7.3.2.

    @ivar stateUpdateFunction: Non-linear function to propagate state into the
        next timestep. The function is called with two arguments, the current
        Nx1 state vector and the current Mx1 control input vector. It should
        return a new Nx1 state vector. N and M are the numbers of states and
        control inputs.
    @ivar measurementFunction: Non-linear function relating state variables
        to expected measurements. The function is called with the current Nx1
        state vector as an argument, and should return an Mx1 measurement
        vector. N and M are the numbers of states and measurement outputs.
    @ivar state: Current state estimate.
    @ivar stateCovariance: State covariance matrix.
    @ivar processCovariance: Process noise covariance matrix.
    @ivar measurementCovariance: Measurement noise covariance matrix.
    """

    def __init__(self, stateUpdateFunction, measurementFunction, state,
            stateCovariance, processCovariance, measurementCovariance):
        """
        Construct unscented Kalman filter.

        @param stateUpdateFunction: Non-linear function to propagate state into
            the next timestep. The function is called with two arguments, the
            current Nx1 state vector and the current Mx1 control input vector.
            It should return a new Nx1 state vector. N and M are the numbers of
            states and control inputs.
        @param measurementFunction: Non-linear function relating state
            variables to expected measurements. The function is called with the
            current Nx1 state vector as an argument, and should return an Mx1
            measurement vector. N and M are the numbers of states and measurement
            outputs.
        @param state: Initial state estimate.
        @param stateCovariance: Initial state covariance matrix.
        @param processCovariance: Process noise covariance matrix.
        @param measurementCovariance: Measurement noise covariance matrix.
        """

        self.stateUpdateFunction = stateUpdateFunction
        self.measurementFunction = measurementFunction

        self.state = state
        self.stateCovariance = stateCovariance
        self.processCovariance = processCovariance
        self.measurementCovariance = measurementCovariance

    @property
    def stateUpdateFunction(self):
        return self._stateUpdateUT._function

    @stateUpdateFunction.setter
    def stateUpdateFunction(self, function):
        self._stateUpdateUT = UnscentedTransform(function)

    @property
    def measurementFunction(self):
        return self._measurementUT._function

    @measurementFunction.setter
    def measurementFunction(self, function):
        self._measurementUT = UnscentedTransform(function)

    def predict(self, control=None):
        state, stateCovariance = self._stateUpdateUT(self._x, self._P, *[control])
        self._x = state
        self._P = stateCovariance + self._Q
        assert self._x.shape == (self._states,1)
        assert np.all(np.isfinite(self._x))
        assert self._P.shape == (self._states, self._states)
        assert np.all(np.isfinite(self._P))

    def innovation(self, measurement):
        measurement = np.asanyarray(measurement).reshape(-1,1)
        predictedMeasurement, predictionCovariance, \
            stateSigmas, measurementSigmas, weights = \
                self._measurementUT(self._x, self._P, returnSigmas=True)
        x,y = self._x, predictedMeasurement
        sigmaPoints = zip(weights, stateSigmas, measurementSigmas)
        innovation = measurement - predictedMeasurement
        innovationCovariance = self._R + predictionCovariance
        crossCovariance = np.sum((w*(X-x)*(Y-y).T for (w,X,Y) in sigmaPoints), axis=0)
        return innovation, innovationCovariance, crossCovariance

    def correct(self, measurement):
        measurement = np.asanyarray(measurement).reshape(-1,1)
        innovation, innovationCovariance, crossCovariance = \
                self.innovation(measurement)
        self._K = np.dot(crossCovariance, inv(innovationCovariance))
        self._x = self._x + self._K * innovation
        assert self._x.shape == (self._states,1)
        assert np.all(np.isfinite(self._x))
        self._P = self._P - self._K * innovationCovariance * self._K.T
        assert self._P.shape == (self._states, self._states)
        assert np.all(np.isfinite(self._P))
