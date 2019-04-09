"""
Spline fitting.
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
from abc import ABC, abstractmethod
from scipy.interpolate import splrep, splev
import math


class Spline(ABC):
    """
    Base class for splines.
    """
    class InsufficientPointsError(ValueError):
        """
        Exception to be raised if insufficient points to form a valid spline.
        """
        pass

    @property
    @abstractmethod
    def validFrom(self):
        """
        The lowest x value at which the spline is defined.
        """
        pass

    @property
    @abstractmethod
    def validTo(self):
        """
        The highest x value atwhich the spline is defined.
        """
        pass

    @abstractmethod
    def __call__(self, x):
        """
        Evaluate the spline output at value(s) x.
        """
        pass


class UnivariateSpline(Spline):
    """
    Model of a function of a single variable using spline fitting.
    """
    def __init__(self, x, y, order=3, stddev=0):
        """
        Construct spline.

        @param x: Length N, monotonically increasing sequence of input values.
        @param y: Length N sequence of output values.
        @param order: Order of spline to fit, 1 <= k <= 5. Odd orders are
            recommended.
        @param stddev: Standard deviation in output values. If non-zero, the
            resulting spline will be smoothed within this constraint.
        """
        splineKwArgs = {'k':order}
        if stddev != 0:
            w = np.zeros_like(x)
            w.fill(1/stddev)
            splineKwArgs['w']=w
        else:
            splineKwArgs['s']=0

        if len(x) <= order:
            raise Spline.InsufficientPointsError("%d points insufficient for order %d spline" % (len(x), order))

        self._tck = splrep(x,y,**splineKwArgs)

        self._validFrom = x[self._boundaryPoints]
        self._validTo = x[-(self._boundaryPoints + 1)]

    @property
    def _boundaryPoints(self):
        return int(math.ceil(self._tck[2]/2))

    @property
    def validFrom(self):
        return self._validFrom

    @property
    def validTo(self):
        return self._validTo

    def __call__(self,x,n=0):
        """
        Evaluate the n-th derivative of the function at input value(s) x.
        """
        return splev(x, self._tck, n)


class PartialInputSpline(Spline):
    """
    Piecewise spline that allows for undefined values in the output domain.

    Splines are generated for each region in which there is valid data.

    This class implements 1D to 1D fitting using L{UnivariateSpline}. It can
    be subclassed to accomodate other output types by overriding the
    _splineClass attribute, and the _validity and _output methods.
    """
    _splineClass = UnivariateSpline

    def __init__(self, x, y, **kwargs):
        """
        Construct partial input spline.

        @param x: Values in the input domain of the spline.
        @param y: Values in the output domain of the spline.
        @param kwargs: Additional arguments to pass to spline constructors.
        """

        # Find indices and x values at which the array becomes (in)valid.
        valid = np.array([False] + list(self._validity(y)) + [False])
        lastValid = np.roll(valid, 1)
        nextValid = np.roll(valid, -1)
        changed = valid ^ lastValid
        becameValid = valid & changed
        starting = becameValid & nextValid
        ending = lastValid & valid & ~nextValid
        starts = list(starting.nonzero()[0] - 1)
        ends = list(ending.nonzero()[0] - 1)

        self.splines = []
        for start, end in zip(starts, ends):
            try:
                xvals = x[start:end+1]
                if len(np.shape(y)) == 1:
                    yvals = y[start:end+1]
                else:
                    yvals = y[:,start:end+1]
                self.splines.append(self._splineClass(xvals,yvals, **kwargs))
            except Spline.InsufficientPointsError:
                starts.remove(start)
                ends.remove(end)

        if len(self.splines) == 0:
            raise Spline.InsufficientPointsError("No valid regions long enough to create a spline")

        xstarts = [s.validFrom for s in self.splines]
        xends = [s.validTo for s in self.splines]
        self.validRegions = list(zip(xstarts, xends))
        self.regions = list(zip(xstarts, xends, self.splines))

    def _validity(self, y):
        """
        Return a boolean array indicating where y is valid.

        This method should be overridden by subclasses using other data types.
        """
        return ~np.isnan(y)

    def _output(self, x, conditions, results, undefined):
        """
        Create combined output array from partial results and default values.

        This method should be overridden by subclasses using other data types.

        @param x: Input values, length N L{np.ndarray}.
        @param conditions: List of length N boolean arrays identifying the
            indices of x for which each set of results should be applied.
        @param results: List of spline results, types as appropriate, with
            lengths corresponding to the number of True values in each
            condition array.
        @param undefined: Boolean array of length N indicating indices of x
            which should be filled with default values.

        @return: Output of length N, type as appropriate.
        """
        out = np.empty_like(np.atleast_1d(x))

        for condition, result in zip(conditions, results):
            out[condition] = result

        out[undefined] = np.nan

        if np.isscalar(x):
            return out[0]
        else:
            return out

    def __call__(self, x, *args, **kwargs):
        X = np.atleast_1d(x)
        conditions = np.atleast_2d([(X>=start) & (X<=end)
            for start,end in self.validRegions])
        undefined = ~np.any(conditions, axis=0)
        if np.all(undefined):
            results = []
            conditions = []
        else:
            results, conditions = zip(*[(s(X[c], *args, **kwargs), c)
                for c,s in zip(conditions, self.splines) if X[c].size > 0])
        return self._output(x, conditions, results, undefined)

    @property
    def validFrom(self):
        return self.splines[0].validFrom

    @property
    def validTo(self):
        return self.splines[-1].validTo
