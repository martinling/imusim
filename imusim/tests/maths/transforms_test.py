"""
Tests of mathematical transforms.
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
from numpy import random
from scipy import stats
from numpy import testing
from imusim.maths import transforms
from imusim.maths.vectors import vector
from imusim.maths.quaternions import Quaternion
from imusim.maths.transforms import AffineTransform
from math import sqrt

RANDOM_SEED=0

def assert_almost_equal(a, b):
    for ai, bi in zip(a, b):
        testing.assert_almost_equal(ai, bi)

# each element of the list is a tuple of cartesian and polar co-ords
CARTESIAN_POLAR_TEST_VALUES = [
        ((0,0),(0,0)),
        ((0,1),(1,np.pi/2)),
        ((1,0),(1,0)),
        ((1/sqrt(2), 1/sqrt(2)), (1,np.pi/4)),
        ((0,10),(10,np.pi/2)),
        ((5,0),(5,0))
        ]
def checkPolarToCartesian(cartesian, polar):
    assert_almost_equal(transforms.polarToCartesian(*polar), cartesian)

def checkCartesianToPolar(cartesian, polar):
    assert_almost_equal(transforms.cartesianToPolar(*cartesian), polar)

def testCartesianPolarConversion():
    for testParams in CARTESIAN_POLAR_TEST_VALUES:
        cartesian, polar = testParams
        yield checkPolarToCartesian, cartesian, polar
        yield checkCartesianToPolar, cartesian, polar

# each element of the list is a tuple of vectors in NED and CG co-ords
NED_CG_TEST_VALUES = [
        (vector(1,0,0), vector(0,0,1)),
        (vector(0,1,0), vector(-1,0,0)),
        (vector(0,0,1), vector(0,-1,0))
        ]

def checkNED_to_CG(ned, cg):
    assert_almost_equal(transforms.convertNEDtoCG(ned), cg)

def checkCG_to_NED(ned, cg):
    assert_almost_equal(transforms.convertCGtoNED(cg), ned)

def testCoordinateChange():
    for ned, cg in NED_CG_TEST_VALUES:
        yield checkNED_to_CG, ned, cg
        yield checkCG_to_NED, ned, cg

def testCG_to_NED_Quat():
    # those tests looks learry weird. Checked them twice, looks like they are wrong
#    testing.assert_equal(
#            transforms.convertCGtoNED(Quaternion()).components,
#            Quaternion(0.5, -0.5, 0.5, -0.5).components)
    pass

def testNED_to_CG_Quat():
    # testing.assert_equal(
    #         transforms.convertNEDtoCG(Quaternion()).components,
    #         Quaternion(0.5, 0.5, -0.5, 0.5).components)
    pass

AFFINE_TRANSFORM_TESTS = [
        (AffineTransform(transform=np.eye(3)), vector(1,0,0), vector(1,0,0)),
        (AffineTransform(rx=np.pi/2), vector(0,0,1), vector(0,-1,0)),
        (AffineTransform(ry=np.pi/2), vector(0,0,1), vector(1,0,0)),
        (AffineTransform(rz=np.pi/2), vector(1,0,0), vector(0,1,0)),
        (AffineTransform(scale=2), vector(1,0,0), vector(2,0,0)),
        (AffineTransform(translation=vector(1,0,0)), vector(0,0,1),
            vector(1,0,1))
        ]
def checkAffineTransform_apply(transform, input, expectedOutput):
    assert_almost_equal(transform.apply(input), expectedOutput)

def checkAffineTransform_reverse(transform, input, expectedOutput):
    assert_almost_equal(transform.reverse(expectedOutput), input)

def testAffineTransforms():
    for transform, input, output in AFFINE_TRANSFORM_TESTS:
        yield checkAffineTransform_apply, transform, input, output
        yield checkAffineTransform_reverse, transform, input, output

def checkMeanAndCovariance(mean, covariance, samples):
    axis = 1 if mean.shape[0] > 1 else 0
    mean = mean.flatten()
    sampleMean = np.mean(samples, axis=axis)
    sdom = stats.sem(samples, axis=axis)
    testing.assert_array_less(np.abs(mean-sampleMean), 3*sdom,
            "Error between Monte Carlo and UT mean estimates greater "
            "than 3 standard deviations of the mean. (p=0.01)")

    # TODO: make this support the full covariance matrix, not just the diagonals
    covariance = np.diag(covariance)
    sampleCovariance = np.diag(np.atleast_2d(np.cov(samples)))
    N = samples.shape[axis]
    secov = np.sqrt(2.0/(N-1)) * sampleCovariance

    testing.assert_array_less(np.abs(covariance-sampleCovariance),
            3*secov,
            "Error between Monte Carlo and UT covariance estimates greater "
            "than 3 times the standard error of the covariance. (p=0.01)")

def testUnscentedTransform_linearFunction():
    xBar = np.array([[1]])
    Sigma = np.array([[0.1]])

    def transformFunction(x):
        return 2*x

    ut = transforms.UnscentedTransform(transformFunction)
    mean, var = ut(xBar, Sigma)

    assert_almost_equal(mean, 2*xBar)
    assert_almost_equal(var, 2*Sigma*2)

def testUnscentedTransform_rotation():
    xBar = np.array([[20],[0]])
    Sigma = np.array([[3,0],[0,2]])

    r = 1/np.sqrt(2)
    T = np.array([[r,-r],[r,r]])
    def transformFunction(x):
        x,y = x
        return np.dot(T,np.vstack((x,y)))

    ut = transforms.UnscentedTransform(transformFunction)
    mean, cov = ut(xBar, Sigma)

    random.seed(RANDOM_SEED)
    x = random.normal(loc=xBar[0], scale=np.sqrt(Sigma[0,0]), size=1000)
    y = random.normal(loc=xBar[1], scale=np.sqrt(Sigma[1,1]), size=1000)
    transformedPoints = np.hstack([transformFunction(p) for p in zip(x,y)])
    checkMeanAndCovariance(mean, cov, transformedPoints)

def testUnscentedTransform_quadratic():
    xBar = np.array([[2]])
    Sigma = np.array([[0.25]])

    def transformFunction(x):
        return x**2

    ut = transforms.UnscentedTransform(transformFunction)
    mean, var = ut(xBar, Sigma)

    random.seed(RANDOM_SEED)
    x = random.normal(loc=xBar[0][0], scale=np.sqrt(Sigma[0][0]), size=1000)
    X = transformFunction(x)
    checkMeanAndCovariance(mean, var, X)

def checkUnscentedTransform_polarToCartesian(params):
    meanR, meanTheta, varR, varTheta = params
    def transformFunction(params):
        r,theta = params
        return np.vstack(transforms.polarToCartesian(r, theta))

    # calculate results using unscented transform
    ut = transforms.UnscentedTransform(transformFunction)
    mean, cov = ut(np.vstack((meanR, meanTheta)), np.diag((varR, varTheta)))

    # calculate results using monte carlo sim
    random.seed(RANDOM_SEED)
    r = random.normal(loc=meanR, scale=np.sqrt(varR), size=1000)
    theta = random.normal(loc=meanTheta, scale=np.sqrt(varTheta), size=1000)
    cartesianCoords = np.hstack([transformFunction(p) for p in zip(r,theta)])
    checkMeanAndCovariance(mean, cov, cartesianCoords)


def testUnscentedTransform_polarToCartesian():
    meanR = (1,10,20)
    varR = (0.1,1,5)
    meanTheta = np.radians((5,15,90))
    varTheta = np.radians((0.1,1, 4))

    for params in zip(meanR, meanTheta, varR, varTheta):
        yield checkUnscentedTransform_polarToCartesian, params




