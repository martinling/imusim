"""
Test for whether ideal Kalman filter process covariances give optimal results.
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

from imusim.maths.kalman import KalmanFilter, UnscentedKalmanFilter
from scipy.signal import lfilter
from scipy.optimize import fmin
import numpy as np

np.random.seed(42)

T = 50
dt = 0.1
t = np.arange(0, T, dt)
c = 0.9
q = 0.01
r = 0.01

def kf(q, states):

	stateUpdate = np.diag([c] * states)
	control = None
	measurement = np.diag([1] * states)
	initialState = np.array([[0]] * states)
	initialCovariance = np.diag([0.5] * states)
	measurementCovariance = np.diag([r] * states)

	return KalmanFilter(
		stateUpdate,
		control,
		measurement,
		initialState,
		initialCovariance,
		np.diag(list(q) * states),
		measurementCovariance)

def ukf(q, states):

	stateUpdateFunc = lambda state, control: state * c
	measurementFunc = lambda state: state
	initialState = np.array([[0]] * states)
	initialCovariance = np.diag([0.5] * states)
	measurementCovariance = np.diag([r] * states)

	return UnscentedKalmanFilter(
		stateUpdateFunc,
		measurementFunc,
		initialState,
		initialCovariance,
		np.diag(list(q) * states),
		measurementCovariance)

def checkOptimalQ(states):
    x = np.array([lfilter([1.0], [1.0, -c], np.random.normal(scale=np.sqrt(q), size=t.shape)) for s in range(states)])
    n = np.random.normal(scale=np.sqrt(r), size=x.shape)
    y = x + n

    for filter in [kf, ukf]:

        def filterError(q, filter):

            filter = filter(q, states)

            z = np.empty_like(x)

            for i in range(len(t)):
                filter.update(y[:,i].reshape(states,1))
                z[:,i] = np.asarray(filter._x)[:,0]

            return np.var(z - x)

        qopt = fmin(filterError, q, (filter,), disp=False)[0]
        np.testing.assert_approx_equal(qopt, q, 1)

def testOptimalQ():
    for states in range(1,4):
        yield checkOptimalQ, states
