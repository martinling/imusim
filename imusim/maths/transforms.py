"""
Mathematical transforms.
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
from imusim.maths.quaternions import Quaternion, QuaternionArray
from imusim.maths.vectors import vector
from imusim.maths.matrices import matrixFromEuler
import collections
import functools
from numpy import linalg

class AffineTransform(object):
    """
    Affine transform composed of a linear transform matrix and translation.
    """
    def __init__(self, transform=None, rz=0, ry=0, rx=0, scale=1,
            translation=vector(0,0,0)):
        """
        Construct an affine transform.

        A generic transformation matrix can be supplied, or a transform
        composed of rotation and scaling can be built.

        Generating rotations using supplied Euler angles is performed from the
        perspective of a fixed camera viewing a rotating object, i.e. points
        rotate within a fixed reference frame. Rotations are applied in
        aerospace (ZYX) order.

        @param transform: 3x3 linear transformation matrix.
        @param rx: X rotation angle in radians.
        @param ry: Y rotation angle in radians.
        @param rz: Z rotation angle in radians.
        @param scale: Scaling factor
        @param translation: 3x1 translation vector.
        """
        if transform is not None:
            self._transform = transform
        else:
            self._transform = scale * matrixFromEuler((rz,ry,rx),'zyx',False)

        self._inverseTransform = np.linalg.inv(self._transform)
        self._translation = translation

    def apply(self,v):
        """
        Apply transform to an array of column vectors.
        """
        return np.tensordot(self._transform, v, axes=([1],[0])) \
                + self._translation

    def reverse(self,v):
        """
        Apply the inverse of this transform to an array of column vectors.
        """
        return np.tensordot(self._inverseTransform,
                v-self._translation,axes=([1],[0]))

class UnscentedTransform(object):
    """
    Implementation of the Unscented Transform of Julier and Uhlmann.

    The unscented transform propagates mean and covariance information
    through a non-linear function. An instance of this class is callable
    to apply the function to a vector of mean values and a covariance
    matrix, returning an output mean vector and covariance matrix.
    """

    _sum = functools.partial(np.sum, axis=0)

    def __init__(self, function):
        """
        Construct unscented transform.

        @param function: The function to apply.
        """
        assert isinstance(function, collections.Callable)
        self._function = function

    def __call__(self, mean, covariance, *args, **kwargs):
        """
        Propagate mean and covariance values through the transform function.

        Additional arguments are passed through to the function.

        @param mean: Nx1 column vector of mean input values.
        @param covariance: NxN covariance of input values.
        @param returnSigmas: If True, return sigma points and weights.

        @return: Mx1 column vector of mean output values, MxM covariance
            of output values. If returnSigmas is True, additionally a list of
            input sigma point column vectors, a list of output sigma point
            column vectors, and a list of sigma point weights.
        """
        returnSigmas = kwargs.pop('returnSigmas', False)
        inputSigmaPoints, weights = UnscentedTransform.sigmaPoints(mean, covariance)
        outputSigmaPoints = [self._function(p, *args, **kwargs) for p in inputSigmaPoints]
        weightedOutputPoints = list(zip(weights, outputSigmaPoints))
        outputMean = self._sum(W*y for W,y in weightedOutputPoints)
        outputCovariance = self._sum(W*(y-outputMean)*(y-outputMean).T for W,y in weightedOutputPoints)
        if returnSigmas:
            return outputMean, outputCovariance, inputSigmaPoints, outputSigmaPoints, weights
        else:
            return outputMean, outputCovariance

    @staticmethod
    def sigmaPoints(mean, covariance):
        """
        Generate sigma points around the given mean values based on covariance.

        Implements the symmetric sigma point set of Julier and Uhlmann 2004,
        'Unscented Filtering and Nonlinear Estimation', equation 12. Weighting
        applied to first sigma point is 1/3 based on the assumption of Gaussian
        distributions.

        @return: List of sigma point column vectors, list of weights.
        """
        N = len(mean)
        mean = np.reshape(mean, (N,1))
        assert covariance.shape == (N,N)

        sigmaPoints = [mean] * (2*N + 1)
        w0 = 1/3 # based on assumption of Gaussian distributions.

        cholesky = linalg.cholesky((N/(1-w0)) * covariance)
        # cholesky returns A s.t. A*A.T = P so we use the columns of A
        columns = np.hsplit(cholesky, N)
        for i, column in enumerate(columns):
            sigmaPoints[i+1] = mean + column
            sigmaPoints[i+1+N] = mean - column
        weights = [w0] + [(1-w0)/(2*N)] * (2 * N)
        return sigmaPoints, weights


def cartesianToPolar(x,y):
    """
    Convert Cartesian co-ordinates (x,y) to polar co-ordinates (r,theta).
    """
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)

    return r,theta

def polarToCartesian(r,theta):
    """
    Convert polar co-ordinates (r,theta) to Cartesian (x,y).
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x,y

# Quaternion for converting between CG and NED co-ordinate frames
_CGtoNED = Quaternion(0.5,0.5,-0.5,0.5)
_NEDtoCG = _CGtoNED.conjugate

def convertCGtoNED(param):
    """
    Convert parameters specified in Computer Graphics co-ordinates to
    North East Down co-ordinates.

    @param param: A L{Quaternion} or 3x1 vector in CG frame.

    @return:  The corresponding L{Quaternion} or vector in NED frame
    """
    if isinstance(param, (Quaternion,QuaternionArray)):
        return _CGtoNED.conjugate * param * _CGtoNED
    else:
        return _CGtoNED.rotateFrame(param)

def convertNEDtoCG(param):
    """
    Convert parameters specified in North East Down co-ordinates to
    Computer Graphics co-ordinates.

    @param param: A L{Quaternion} or 3x1 vector in NED frame.

    @return:  The corresponding L{Quaternion} or vector in CG frame
    """
    if isinstance(param, (Quaternion,QuaternionArray)):
        return _NEDtoCG.conjugate * param * _NEDtoCG
    else:
        return _NEDtoCG.rotateFrame(param)
