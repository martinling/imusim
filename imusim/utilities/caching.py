"""
Decorators to enable caching of function results.
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

from functools import wraps
from collections import defaultdict

class CacheEntry(object):
    __slots__ = ('time', 'result')

    def __init__(self):
        self.time = -1
        self.result = None

class CacheLastValue(object):
    """
    Provides a decorator for methods of the form f(t)
    to cache the last value

    @param tolerance: tolerance to use when checking that two times are equal.
    """

    def __init__(self, tolerance=1e-6):
        self._tolerance = tolerance
        self._cache = defaultdict(CacheEntry)

    def __call__(self, method):
        tolerance = self._tolerance
        closure_abs = abs
        cache = self._cache

        @wraps(method)
        def checkcache(obj, t):
            cacheEntry = cache[obj]
            try:
                if closure_abs(t - cacheEntry.time) < tolerance:
                    return cacheEntry.result
                cacheEntry.time = t
                cacheEntry.result = method(obj, t)
                return cacheEntry.result
            except:
                return method(obj, t)
        return checkcache


