"""
Tests for integration methods.
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
from imusim.maths import integrators
from imusim.testing.vectors import assert_vectors_correlated
from imusim.testing.inspection import getImplementations

x = np.linspace(0,10,100)
dt = x[1]-x[0]

functions = ( lambda x: x,
        lambda x: x**2,
        lambda x: np.cos(x),
        lambda x: np.e**x )
integrals = ( lambda x: 1/2 * x**2,
        lambda x: 1/3 * x**3,
        lambda x: np.sin(x),
        lambda x: np.e**x)
doubleIntegrals = (lambda x: x**3 / 6,
        lambda x: x**4 / 12,
        lambda x: -np.cos(x),
        lambda x: np.e**x)

def checkIntegral(method, function, integral):
    estimated_integral = np.empty_like(x)
    integrator = method(integral(0))
    for i,s in enumerate(x[1:]):
        estimated_integral[i+1] = integrator(function(s), dt)
    assert_vectors_correlated(estimated_integral[1:], integral(x[1:]))

def checkDoubleIntegral(method, function, integral, doubleIntegral):
    estimated_integral = np.empty_like(x)
    integrator = integrators.DoubleIntegrator(doubleIntegral(0),integral(0),
            method)
    for i,s in enumerate(x[1:]):
        estimated_integral[i+1] = integrator(function(s), dt)
    assert_vectors_correlated(estimated_integral[1:], doubleIntegral(x[1:]))

def testIntegrators():
    for method in getImplementations(integrators, integrators.Integrator):
        for function, integral, doubleIntegral in zip(functions, integrals,
                doubleIntegrals):
            yield checkIntegral, method, function, integral
            yield (checkDoubleIntegral, method, function, integral,
                doubleIntegral)

