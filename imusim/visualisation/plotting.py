"""
Support for plotting IMUSim data types.
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
import pylab
import numpy as np
import numbers
from imusim.utilities.time_series import TimeSeries
from imusim.maths.quaternions import Quaternion, QuaternionArray

def plotTimeSeries(timeSeries, *args, **kwargs):
    """
    Plot data contained in a time series.

    Additional *args and **kwargs are passed on to pylab.plot.
    """
    x = timeSeries.timestamps
    y = timeSeries.values
    var = timeSeries.variances
    label = kwargs.get('label', '')
    if issubclass(timeSeries.dtype, np.ndarray):
        if y.shape[0] == 3:
            # We will assume these are vectors
            for i, axis in enumerate(('x', 'y', 'z')):
                kwargs['label'] = ' '.join((label, axis))
                plotWithVariance(x, y[i], var[:,i,i], *args, **kwargs)
        else:
            if y.ndim > 1:
                for i in range(len(y)):
                    plotWithVariance(x, y[i], var[:,i,i], *args, **kwargs)
            else:
                plotWithVariance(x, y, var, *args, **kwargs)

    elif issubclass(timeSeries.dtype, Quaternion):
        for i, component, in enumerate(('w', 'x', 'y', 'z')):
            kwargs['label'] = ' '.join((label, component))
            plotWithVariance(x, y[:,i], var[:,i,i], *args, **kwargs)

    elif issubclass(timeSeries.dtype, numbers.Real):
        plotWithVariance(x, y, var, *args, **kwargs)

def plotWithVariance(x, y, variance, *args, **kwargs):
    """
    Plot data with variance indicated by shading within one sigma.
    """
    line = pylab.plot(x, y.flatten(), *args, **kwargs)[0]
    sigma = np.sqrt(variance)
    pylab.fill_between(x, y - sigma, y + sigma, color=line.get_color(), alpha=0.5)

def plot(x, y=None, *args, **kwargs):
    """
    Smart version of pylab.plot(...) that understands types used in IMUsim.

    The supported y value types are:

        - L{TimeSeries} data
        - Lists or tuples of Quaternions
        - Lists or tuples of column vectors
        - Arrays of column vectors
    """

    singleArg = (y is None)

    if singleArg:
        if isinstance(x, TimeSeries):
            plotTimeSeries(x, *args, **kwargs)
        else:
            y = x

    if isinstance(y, (list,tuple)):
        if isinstance(y[0], Quaternion):
            y = QuaternionArray(y).array
        elif isinstance(y[0], np.ndarray):
            y = np.hstack(y).T
    elif isinstance(y, QuaternionArray):
        y = y.array
    elif isinstance(y, np.ndarray):
        if y.ndim == 2:
            y = y.T

    # if singleArg:
    #     pylab.plot(y, *args, **kwargs)
    # else:
    #     pylab.plot(x, y, *args, **kwargs)
