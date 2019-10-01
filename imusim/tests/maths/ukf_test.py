"""
Tests for unscented Kalman filter implementation.
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
from imusim.maths.kalman import UnscentedKalmanFilter

def testNonLinearMeasurement():
    """
    Tests filter estimating constants with a non-linear measurement function
    """
    SAMPLES = 1000

    stateUpdate = lambda state, control: state

    def measurementFunc(state):
        state = np.asarray(state)
        newState = np.empty_like(state)
        newState[0] = state[0]**2
        newState[1] = state[1]**3
        return newState

    initialState = np.matrix('0.1;0.1')
    initialCovariance = np.matrix('1,0;0,1')
    processCovariance = np.matrix('0.001,0;0,0.001')
    measurementCovariance = np.matrix('0.01,0;0,0.01')

    ukf = UnscentedKalmanFilter(stateUpdate, measurementFunc, initialState,
            initialCovariance, processCovariance, measurementCovariance)

    trueStates = np.empty((2,SAMPLES))
    trueStates.fill(2)
    measurements = measurementFunc(trueStates) + np.random.normal(loc=0,
            scale=0.1,
            size=(2,SAMPLES))
    estimates = np.empty_like(trueStates)
    estimateCovariances = np.empty((SAMPLES,2,2))

    for i in range(SAMPLES):
        ukf.update(measurements[:,i:i+1])
        estimates[:,i:i+1] = ukf.state
        estimateCovariances[i] = ukf.stateCovariance

    # Allow ten samples for estimates to start settling
    errors = trueStates[:,10:]-estimates[:,10:]
    assert np.mean(errors[0]**2) < 0.01
    assert np.mean(errors[1]**2) < 0.01
    assert np.var(errors[0]) < measurementCovariance[0,0]
    assert np.var(errors[1]) < measurementCovariance[1,1]

def testNonLinearStateUpdate():
    """
    Test filter with a non-linear state update
    """

    SAMPLES = 1000

    def stateUpdate(state, control):
        state = np.asarray(state)
        newState = np.empty_like(state)
        newState[0] = state[0] + 0.001*state[1]**2
        newState[1] = state[1] + 0.01
        return newState

    measurementFunc = lambda state: state

    initialState = np.matrix('0;0')
    initialCovariance = np.matrix('1,0.5;0.5,1')
    processCovariance = np.matrix('0.01,0.01;0.01,0.01')
    measurementCovariance = np.matrix('0.25,0;0,0.25')

    ukf = UnscentedKalmanFilter(stateUpdate, measurementFunc, initialState,
            initialCovariance, processCovariance, measurementCovariance)

    trueStates = np.empty((2,SAMPLES))
    trueStates[:,0] = np.array([1,2])
    measurements = np.empty_like(trueStates)
    estimates = np.empty_like(trueStates)
    estimateCovariances = np.empty((SAMPLES,2,2))

    for i in range(1,SAMPLES):
        trueStates[:,i] = stateUpdate(trueStates[:,i-1],0)
        measurements[:,i] = trueStates[:,i] + np.random.normal(loc=0, scale=0.5,
                size=2)

        ukf.update(measurements[:,i:i+1])

        estimates[:,i:i+1] = ukf.state
        estimateCovariances[i] = ukf.stateCovariance

    # Allow ten samples for estimates to start settling
    errors = trueStates[:,10:]-estimates[:,10:]
    assert np.mean(errors[0]**2) < 0.1
    assert np.mean(errors[1]**2) < 0.1
    assert np.var(errors[0]) < measurementCovariance[0,0]
    assert np.var(errors[1]) < measurementCovariance[1,1]

