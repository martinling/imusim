"""
Functions for working with column vectors.
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
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)

__all__ = ['vector']

def vector(*components):
    return np.array([components],dtype=np.float64).T

def norm(np.ndarray[np.float64_t, ndim=2] vectors):
    """
    Compute the l2norm of an array of column vectors
    """
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef np.ndarray[np.float64_t, ndim=1] norms = np.empty(vectors.shape[1])
    for i in range(vectors.shape[1]):
        norms[i] = 0
        for j in range(vectors.shape[0]):
            norms[i] += vectors[j,i] * vectors[j,i]
        norms[i] = sqrt(norms[i])
    return norms

def cross(np.ndarray[np.float64_t, ndim=2] v1,
        np.ndarray[np.float64_t, ndim=2] v2):
    """
    Compute the cross products of two arrays of column vectors.

    Each array of vectors must have the same shape
    """
    cdef Py_ssize_t i
    cdef np.ndarray[np.float64_t, ndim=2] crossProduct = np.empty((3,v1.shape[1]))

    assert v1.shape[0] == v2.shape[0] == 3 and v1.shape[1] == v2.shape[1], \
            "Vectors must both have shape (3,%d). (v1.shape=(%d, %d), v2." \
            "shape=(%d, %d)" %(v1.shape[1], v1.shape[0], v1.shape[1],
                    v2.shape[0], v2.shape[1])

    for i in range(v1.shape[1]):
        crossProduct[0,i] = v1[1,i] * v2[2,i] - v1[2,i] * v2[1,i]
        crossProduct[1,i] = v1[2,i] * v2[0,i] - v1[0,i] * v2[2,i]
        crossProduct[2,i] = v1[0,i] * v2[1,i] - v1[1,i] * v2[0,i]
    return crossProduct

def dot(np.ndarray[np.float64_t, ndim=2] v1,
        np.ndarray[np.float64_t, ndim=2] v2):
    """
    Compute the dot products of two arrays of column vectors.

    Each array must have the same shape
    """
    cdef Py_ssize_t i
    cdef np.ndarray[np.float64_t, ndim=1] dotProduct = np.empty(v1.shape[1])

    assert v1.shape[0] == v2.shape[0] and v1.shape[1] == v2.shape[1], \
            "Vectors must both have shape (%d,%d). (v1.shape=(%d, %d), v2." \
            "shape=(%d, %d)" %(v1.shape[0], v1.shape[1], v1.shape[0], v1.shape[1],
                    v2.shape[0], v2.shape[1])

    for i in range(v1.shape[1]):
        dotProduct[i] = 0
        for j in range(v1.shape[0]):
            dotProduct[i] += v1[j,i] * v2[j,i]
    return dotProduct

def validity(vectors):
    """
    Return a boolean array indicating which vectors in an array of column
    vectors do not contain NaN values.
    """
    return ~np.any(np.isnan(vectors), axis=0)

def nan(length=1):
    v = np.empty((3,length))
    v[:] = np.nan
    return v
