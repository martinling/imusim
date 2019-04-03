"""
Classes for representing time series data.
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
from collections import namedtuple
import numbers
from copy import copy
from imusim.maths.quaternions import Quaternion, QuaternionArray

class TimeSeries(object):
    """
    A series of values taken at consecutive points in time.

    The values may be scalars, column vectors, or L{Quaternion} objects.

    Each value may have a (co)variance associated with it.

    Timestamps may be irregularly spaced but must be monotonically increasing.
    """

    def __init__(self, timestamps=None, values=None, variances=None):
        """
        Construct time series.

        A time series can be constructed empty, or with initial values if
        passed to the contstructor.

        @param timestamps: Times at which values were taken.
        @param values: Either a length N array of scalars or other data types,
            an MxN column vector array, or a length N L{QuaternionArray}.
        @param variances : Variances, either a length N sequence of variances
            for scalar data, a NxMxM sequence of covariances for MxN column
            vector data, or an Nx4x4 sequence of covariances for quaternion
            data.
        """
        if (timestamps is None) != (values is None):
            raise ValueError("Both or neither of timestamps and values must be provided")
        self._dtype = None
        self._timestampsArray = None
        self._valuesArray = None
        self._variancesArray = None
        if timestamps is not None:
            assert all(np.diff(timestamps) > 0), "Timestamps must be in order"
            if len(np.shape(values)) == 1:
                values = list(np.asarray(values))
            else:
                values = list(np.asarray(np.hsplit(values, values.shape[-1])))
            if variances is not None:
                [self.add(t,v,V) for t,v,V in zip(timestamps, values, variances)]
            else:
                [self.add(t,v) for t,v in zip(timestamps, values)]

    def __getstate__(self):
        if self._hasVariances:
            return (self.timestamps, self.values, self.variances)
        else:
            return (self.timestamps, self.values, None)

    def __setstate__(self, state):
        self.__init__(*state)

    def __len__(self):
        if self._dtype is None:
            return 0
        else:
            return len(self._timestamps)

    def _initTypes(self, value):
        self._dtype = type(value)
        self._dshape = np.shape(value)
        if isinstance(value, np.ndarray):
            value = np.asarray(value)
            self._dims = value.shape[0]
        elif isinstance(value, Quaternion):
            self._dims = 4
        else:
            self._dims = 1

    def _varianceShape(self):
        if self._dims == 1:
            return ()
        else:
            return (self._dims, self._dims)

    @property
    def dtype(self):
        """
        Type of data in this time series.
        """
        return self._dtype

    @property
    def timestamps(self):
        """
        Array of times at which the values were taken.
        """
        if self._timestampsArray is None:
            self._timestampsArray = np.array(self._timestamps)
        return self._timestampsArray

    def _makeValuesArray(self, values):
        if self.dtype is Quaternion:
            return QuaternionArray(values)
        elif issubclass(self._dtype, np.ndarray) and len(self._dshape) == 2 and self._dshape[1] == 1:
            return np.hstack(values)
        else:
            return np.array(values)

    @property
    def values(self):
        """
        Array of values of the time series in chronological order.

        The returned array has shape (N) for scalars, shape (n, N) for
        vectors, or is a L{QuaternionArray} of length N for quaternions,
        where N is the number of values.
        """
        if self._valuesArray is None:
            self._valuesArray = self._makeValuesArray(self._values)
        return self._valuesArray

    @property
    def variances(self):
        """
        The (co)variances of the values in chronological order.

        The returned array has shape (N) if the values are scalar,
        (N, n, n) for vectors, or (N, 4, 4) for quaternions.
        """
        if not self._hasVariances:
            shape = tuple([len(self)] + list(self._varianceShape()))
            variances = np.empty(shape)
            variances[:] = np.nan
            return variances
        if self._variancesArray is None:
            self._variancesArray = np.array(self._variances)
        return self._variancesArray

    def __call__(self, t, returnVariance=False):
        """
        Obtain the closest previous value to time t.

        If returnVariance is True then the value and variance are returned.
        """
        indices = np.searchsorted(self.timestamps, t, 'right') - 1
        if issubclass(self._dtype, np.ndarray) and len(self._dshape) == 2 and self._dshape[1] == 1:
            values = self.values[:,np.atleast_1d(indices)]
        else:
            values = self.values[indices]
        if not returnVariance:
            return values
        else:
            if not self._hasVariances:
                raise ValueError("This time series has no variance data.")
            return values, np.array(self._variances[indices])

    def __iter__(self):
        return zip(self._timestamps, self._values, self._variances)

    def add(self, time, value, variance=None):
        """
        Add a timestamped value to the data.

        @param time: Time at which the value was taken.
        @param value: Value to store. The value must support the __copy__
            operation as used by the copy module.
        @param variance: (Co)variance of the value. The variance must support
            the __copy__ operation as used by the copy module.
        """
        if isinstance(value, np.matrix):
            value = np.asarray(value)
        if np.shape(value) == (1,):
            value = value[0]
        if self._dtype is None:
            self._initTypes(value)
            self._timestamps = [time]
            self._values = [value]
            if variance is not None:
                self._hasVariances = True
                assert np.shape(variance) == self._varianceShape(), \
                        "Variance shape must match series dimensions"
                self._variances = [variance]
            else:
                self._hasVariances = False
        else:
            assert time > self.latestTime, \
                    "Time must be after existing entries"
            self._timestamps.append(copy(time))
            assert type(value) == self._dtype, \
                    "Data type must match series data type"
            assert np.shape(value) == self._dshape, \
                    "Data shape must match series data shape"
            if self._dims > 1 and self._dtype != Quaternion:
                value = np.asarray(value)
            self._values.append(copy(value))
            if self._hasVariances:
                assert np.shape(variance) == self._varianceShape(), \
                        "Variance shape must match series dimensions"
                if self._dims > 1:
                    variance = np.asarray(variance)
                self._variances.append(copy(variance))

        self._latestTime = time
        self._latestValue = value
        self._timestampsArray = None
        self._valuesArray = None
        self._variancesArray = None

    @property
    def earliestTime(self):
        """ The earliest time stamp in this time series. """
        return self._timestamps[0]

    @property
    def earliestValue(self):
        """ The earliest value in this time series. """
        return self._values[0]

    @property
    def latestTime(self):
        """ The latest time stamp in this time series. """
        return self._timestamps[-1]

    @property
    def latestValue(self):
        """ The latest value in this time series. """
        return self._values[-1]
