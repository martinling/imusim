"""
Utilities for testing quaternion values
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
from imusim.maths.quaternions import Quaternion, QuaternionArray

def assertQuaternionAlmostEqual(q1, q2, tol=1e-6, errMsg=""):
    """
    Test whether two quaternions are equal within a given tolerance.
    """
    if isinstance(q1, Quaternion):
        assert(isinstance(q2, Quaternion))
        assert(np.allclose(q1.components, q2.components, atol=tol) or
                np.allclose((-q1).components, q2.components, atol=tol)), \
                "Quaternions not almost equal: %r <> %r" %(q1, q2) +'\n' + \
                errMsg
    elif isinstance(q1, QuaternionArray):
        assert(isinstance(q2, QuaternionArray))
        assert(np.allclose(q1.array, q2.array, atol=tol) or
                np.allclose((-q1).array, q2.array, atol=tol)), \
                "QuaternionArrays not almost equal: %r <> %r" %(q1, q2) + \
                '\n' + errMsg
    else:
        assert False, "Both arguments must be Quaternions or QuaternionArrays." \
                + '\n' + errMsg

def assertMaximumErrorAngle(q1, q2, maximumError, errMsg=""):
    """
    Assert that the error angle between two quaternions is below a maximum.
    """
    if isinstance(q1, Quaternion):
        assert(isinstance(q2, Quaternion))
        error = np.degrees(2*np.arccos((q1*q2.conjugate).w))
        assert error < maximumError, \
                "Angle between quaternions %.2f > %.2f degrees" %(
                    error, maximumError) +'\n' + errMsg
    elif isinstance(q1, QuaternionArray):
        assert(isinstance(q2, QuaternionArray))
        error = np.degrees(2*np.arccos((q1*q2.conjugate).w))
        assert max(error) < maximumError, \
                "Maximum angle between quaternions %.2f > %.2f degrees" %(
                    max(error), maximumError) +'\n' + errMsg
    else:
        assert False, "Both arguments must be Quaternions or QuaternionArrays."\
                +'\n' + errMsg

def assertMaximumRMSErrorAngle(q1, q2, maxRMSE, errMsg=""):
    """
    Assert that max RMSE between 2 L{QuaternionArray}s is less than maxRMSE.

    @param q1: The first L{QuaternionArray}.
    @param q2: The second L{QuaternionArray}.
    @param maxRMSE: The maximum allowable RMS error in degrees.
    @param errMsg: Additional errod detail to display if the assertion fails.
    """
    assert isinstance(q1, QuaternionArray)
    assert isinstance(q2, QuaternionArray)
    errors = np.degrees(2*np.arccos(np.clip(np.abs((q1*q2.conjugate).w),0,1)))
    rmse = np.sqrt(np.mean(errors**2))
    assert rmse < maxRMSE, 'RMSE between quaternions greater than specified' \
            " %.2f>%.2f\n%s"%(rmse, maxRMSE, errMsg)

def assert_quaternions_correlated(actual, desired, targetCorrelation=0.95,
        errMsg=""):
    """
    Assert that correlation between two quaternion arrays is greater than a
    specified target.

    Correlations are calculated as the diagonal of the upper right quadrant
    of the full correlation matrix. Each of the diagonal elements of this
    sub-matrix represents the correlation between an actual vector component
    and its corresponding desired value.

    @param actual: Actual L{QuaternionArray}.
    @param desired: Desired L{QuaternionArray}.
    @param targetCorrelation: Minimum correlation to accept.
    @param errMsg: Additional error detail to display if an assertion fails.
    """
    actual = actual.array.T
    desired = desired.array.T
    correlationMatrix = np.corrcoef(actual, desired)
    s = correlationMatrix.shape[0] // 2
    correlationMatrix = correlationMatrix[s:, :s]
    correlations = np.diag(correlationMatrix)

    assert np.all(correlations > targetCorrelation), "Quaternions are not " \
    "sufficiently correlated. \nCorrelation values: %s\nTarget %.2f"\
        %(correlations, targetCorrelation) + '\n' + errMsg
