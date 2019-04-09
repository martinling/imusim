"""
Iterative integrator implementations.
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

from imusim.utilities.documentation import prepend_method_doc
from abc import ABC, abstractmethod
from copy import copy


class Integrator(ABC):
    """
    Base class for integrators.
    """
    def __init__(self, initialValue):
        """
        Initialise integrator.

        @param initialValue: Initial value for integral.
        """
        self._accumulator = copy(initialValue)

    @abstractmethod
    def __call__(self, sampleValue, dt):
        """
        Update the integral with a new sample value and time step.

        @return: The updated integral value.
        """
        pass


class RectangleRule(Integrator):
    """
    Integration by the rectangle rule.
    """

    def __call__(self, sampleValue, dt):
        self._accumulator += sampleValue * dt
        return copy(self._accumulator)


class TrapeziumRule(Integrator):
    """
    Integration by the trapezium rule.
    """
    def __init__(self, initialValue):
        Integrator.__init__(self, initialValue)
        self._previousValue = None

    def __call__(self, sampleValue, dt):
        if self._previousValue is not None:
            self._accumulator += dt * (self._previousValue + sampleValue) / 2
        else:
            self._accumulator += sampleValue * dt
        self._previousValue = sampleValue
        return copy(self._accumulator)


class DoubleIntegrator(object):
    """
    Iterative integrator that performs double integration.
    """
    @prepend_method_doc(Integrator)
    def __init__(self, initialValue, initialDerivative, method=RectangleRule):
        """
        @param initialDerivative: Initial value for derivative of integral.
        @param method: L{Integrator} class to use.
        """
        self.integrator_1 = method(initialDerivative)
        self.integrator_2 = method(initialValue)

    def __call__(self, sampleValue, dt):
        return self.integrator_2(self.integrator_1(sampleValue, dt), dt)
