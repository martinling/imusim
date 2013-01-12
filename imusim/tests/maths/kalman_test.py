"""
Kalman filter implementation tests.
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
from imusim.maths import kalman
from numpy import testing

t, dt = np.linspace(0, 10, 1000, retstep=True)

stateTransitionMatrices = [
        np.matrix([1]),
        np.matrix([1]),
        np.matrix([[1,dt],[0,1]])
        ]
measurmentMatrices = [
        np.eye(1),
        np.eye(1),
        np.matrix([1,0])
        ]
controlMatrices = [
        np.eye(1),
        np.zeros((1,1)),
        np.matrix([[dt**2/2],[dt]])
        ]
controlFunctions = [
        lambda t: 0.1*t,
        lambda t: 0*t,
        lambda t: np.ones_like(t)
        ]
initialStates = [
        np.zeros((1,1)),
        np.array([[1]]),
        np.array([[0],[0]])
        ]
initialCovariances = [
        np.array([[0.1]]),
        np.array([[0.5]]),
        0.2**2 * np.array([[dt**4/4,dt**3/2],[dt**3/2, dt**2]])
        ]
processCovariances = [
        np.array([[0.01]]),
        np.array([[0.05]]),
        0.2**2 * np.array([[dt**4/4,dt**3/2],[dt**3/2, dt**2]])
        ]
measurementCovariances = [
        np.array([[0.3]]),
        np.array([[0.2]]),
        np.array([[10]])
        ]

def checkKalmanFilter(stateTransition, measurement, control, initialState,
        initialCovariance, processCovariance, measurementCovariance,
        controlFunction):

    states = initialState.shape[0]
    measurements = measurement.shape[0]
    controls = control.shape[1]

    filter = kalman.KalmanFilter(stateTransition, control,
            measurement, initialState, initialCovariance, processCovariance,
            measurementCovariance)

    testing.assert_equal(filter.stateTransitionMatrix, stateTransition)
    testing.assert_equal(filter.controlMatrix, control)
    testing.assert_equal(filter.measurementMatrix, measurement)
    testing.assert_equal(filter.measurementCovariance, measurementCovariance)
    testing.assert_equal(filter.processCovariance, processCovariance)

    signal = np.empty((states, 1000))
    est = np.empty_like(signal)
    meas = np.zeros((measurements, 1000))
    con = controlFunction(t).reshape(controls, 1000)
    signal[:,0:1] = initialState

    processDeviation = np.asmatrix(np.sqrt(processCovariance))
    measurementDeviation = np.asmatrix(np.sqrt(measurementCovariance))

    for i in range(1,1000):
        signal[:,i:i+1] = (stateTransition * signal[:,i-1:i] +
                control * con[:,i] + processDeviation * np.random.normal(
                size=(states, 1)) )
        meas[:,i:i+1] = (measurement * signal[:,i:i+1] + measurementDeviation *
                np.random.normal(size=(measurements, 1)))

        filter.predict(con[:,i])
        filter.correct(meas[:,i])

        est[:,i:i+1] = filter.state

    assert np.all(np.var(est-signal, axis=1) < np.var(meas-signal, axis=1))
    return signal.T, meas.T, est.T

def testKalmanFilter():
    np.random.seed(0)
    for params in zip(stateTransitionMatrices, measurmentMatrices,
            controlMatrices, initialStates, initialCovariances,
            processCovariances, measurementCovariances, controlFunctions):
        yield (checkKalmanFilter, ) + params

def drawFigures():
    import pylab
    for v in testKalmanFilter():
        pylab.figure()
        s,m,e = v[0](*v[1:])
        pylab.plot(t, s, label='signal')
        pylab.plot(t, m, label='meas')
        pylab.plot(t, e, label='est')
        pylab.legend()
    pylab.show()

if __name__ == '__main__':
    drawFigures()



