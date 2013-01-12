"""
Cython wrapper for natural neighbour C implementation.
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

cdef extern from "natural_neighbour/delaunay.h":
    ctypedef struct vertex:
        pass
    ctypedef struct mesh:
        pass
    vertex *initPoints(double *x, double *y, double *z, double *u, double *v, double *w, int n)
    mesh *newMesh()
    void buildMesh(vertex *v, int n, mesh *m)

cdef extern from "natural_neighbour/natural.h":
    void interpolate3_3(double x, double y, double z, double *u, double *v, double *w, mesh *m)

ctypedef np.float64_t dtype_t

cdef class NaturalNeighbourInterpolatorC:

    cdef mesh *m

    def __cinit__(self,
        np.ndarray[dtype_t, ndim=1, mode="c"] x,
        np.ndarray[dtype_t, ndim=1, mode="c"] y,
        np.ndarray[dtype_t, ndim=1, mode="c"] z,
        np.ndarray[dtype_t, ndim=1, mode="c"] u,
        np.ndarray[dtype_t, ndim=1, mode="c"] v,
        np.ndarray[dtype_t, ndim=1, mode="c"] w):

        cdef int n = len(x)
        cdef vertex *points = initPoints(
            <dtype_t *>x.data,
            <dtype_t *>y.data,
            <dtype_t *>z.data,
            <dtype_t *>u.data,
            <dtype_t *>v.data,
            <dtype_t *>w.data,
            n)
        self.m = newMesh()
        buildMesh(points, n, self.m)

    def __call__(self,
        np.ndarray[dtype_t, ndim=1] x,
        np.ndarray[dtype_t, ndim=1] y,
        np.ndarray[dtype_t, ndim=1] z):

        cdef np.ndarray result = np.empty((3,len(x)))

        cdef double U
        cdef double V
        cdef double W
        cdef int i
        for i in range(len(x)):
            interpolate3_3(
                x[i],
                y[i],
                z[i],
                &U,
                &V,
                &W,
                self.m)
            result[0,i] = U
            result[1,i] = V
            result[2,i] = W

        return result
