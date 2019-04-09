"""
Quaternion maths.
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
from scipy import interpolate
from imusim.maths.vector_splines import PartialInputVectorSpline
from imusim.maths import vectors, matrices
import operator
import functools

cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)
    double cos(double x)
    double acos(double x)
    double sin(double x)
    double log(double x)
    double exp(double x)

cdef double EPS = 1e-6

cdef inline void mult_quat_quat(quaternion_t *q, quaternion_t *p,
        quaternion_t *dest):
    cdef double w
    cdef double x
    cdef double y
    cdef double z
    w = q.w*p.w - q.x*p.x - q.y*p.y - q.z*p.z
    x = q.w*p.x + q.x*p.w + q.y*p.z - q.z*p.y
    y = q.w*p.y - q.x*p.z + q.y*p.w + q.z*p.x
    z = q.w*p.z + q.x*p.y - q.y*p.x + q.z*p.w
    dest.w = w
    dest.x = x
    dest.y = y
    dest.z = z

cdef inline void mult_quat_scalar(quaternion_t *q, double scalar,
        quaternion_t *dest):
    cdef double w
    cdef double x
    cdef double y
    cdef double z
    w = q.w*scalar
    x = q.x*scalar
    y = q.y*scalar
    z = q.z*scalar
    dest.w = w
    dest.x = x
    dest.y = y
    dest.z = z

cdef inline void quaternion_add(quaternion_t *q, quaternion_t *p,
        quaternion_t *dest):
    dest.w = q.w + p.w
    dest.x = q.x + p.x
    dest.y = q.y + p.y
    dest.z = q.z + p.z

cdef inline void quaternion_sub(quaternion_t *q, quaternion_t *p,
        quaternion_t *dest):
    dest.w = q.w - p.w
    dest.x = q.x - p.x
    dest.y = q.y - p.y
    dest.z = q.z - p.z

cdef inline void quaternion_log(quaternion_t *q, quaternion_t *dest):
    """
    Obtain the natural logarithm of q.
    """
    cdef double V2 = q.x*q.x + q.y*q.y + q.z*q.z
    cdef double normV = sqrt(V2)
    cdef double m, s, a
    if normV > EPS:
        m = sqrt(q.w**2 + V2)
        a = acos(q.w/m)
        s = a / normV
        dest.w = log(m)
        dest.x = q.x * s
        dest.y = q.y * s
        dest.z = q.z * s
    else:
       dest.w = 0
       dest.x = 0
       dest.y = 0
       dest.z = 0

cdef inline void quaternion_exp(quaternion_t *q, quaternion_t *dest):
    """
    Obtain the exponential of this quaternion.
    """
    cdef double V2 = q.x*q.x + q.y*q.y + q.z*q.z
    cdef double normV = sqrt(V2)
    cdef double s, e
    if normV > EPS:
        s = sin(normV) / normV
        e = exp(q.w)
        dest.w = cos(normV)
        dest.x = q.x * s
        dest.y = q.y * s
        dest.z = q.z * s
        mult_quat_scalar(dest, e, dest)
    else:
        dest.w = exp(q.w)
        dest.x = 0
        dest.y = 0
        dest.z = 0

def QuaternionFromEuler(angles,order='zyx',inDegrees=True):
    """
    Construct a quaternion from an Euler angle sequence.

    @param angles: Sequence of 3 Euler angles.
    @param order: The order to apply the Euler angle sequence.
    @param inDegrees: True to indicate that angles are in degrees (default),
        or False for radians.
    """
    q = Quaternion()
    q.setFromEuler(angles,order,inDegrees)
    return q

def QuaternionFromVectors(x,y,z):
    """
    Create a quaternion that acts as a rotation taking the e1,e2,e3 basis
    vectors to x,y,z.
    """
    q = Quaternion()
    q.setFromVectors(x,y,z)
    return q

def QuaternionFromMatrix(m):
    """
    Create a quaternion from a rotation matrix.
    """
    q = Quaternion()
    q.setFromMatrix(m)
    return q

def QuaternionFromAxisAngle(axis, angle):
    """
    Create a quaternion from an axis and angle.
    """
    return Quaternion(np.cos(angle/2.),
            *(axis.flatten() * np.sin(angle/2.)))

def QuaternionNaN():
    return Quaternion(np.nan, np.nan, np.nan, np.nan)

cdef class Quaternion:
    """
    A quaternion value.
    """
    def __cinit__(self,double w=1,double x=0,double y=0,double z=0):
        self._components.w = w
        self._components.x = x
        self._components.y = y
        self._components.z = z

    # Quaternions are mutable so can't be used for dictionary keys
    __hash__ = None

    property w:
        def __get__(self):
            return self._components.w
        def __set__(self, double value):
            self._components.w = value
    property x:
        def __get__(self):
            return self._components.x
        def __set__(self, double value):
            self._components.x = value
    property y:
        def __get__(self):
            return self._components.y
        def __set__(self, double value):
            self._components.y = value
    property z:
        def __get__(self):
            return self._components.z
        def __set__(self, double value):
            self._components.z = value

    def __richcomp__(Quaternion self, Quaternion other, int op):
        if op == 2:
            return self.w == other.w and \
                self.x == other.x and \
                self.y == other.y and \
                self.z == other.z
        elif op == 3:
            return ~(self.w == other.w and \
                self.x == other.x and \
                self.y == other.y and \
                self.z == other.z)
        else:
            return NotImplemented

    def __repr__(Quaternion self):
        return "Quaternion(%r, %r, %r, %r)"%(self.w,self.x,self.y,self.z)

    property components:
        def __get__(Quaternion self):
            return (self.w,self.x,self.y,self.z)

    def __add__(Quaternion q, Quaternion p):
        return Quaternion(q.w+p.w,q.x+p.x,q.y+p.y,q.z+p.z)

    def __sub__(Quaternion q, Quaternion p):
        return Quaternion(q.w-p.w,q.x-p.x,q.y-p.y,q.z-p.z)

    def __neg__(Quaternion self):
        return Quaternion(-self.w,-self.x,-self.y,-self.z)

    def __mul__(q,p):
        cdef Quaternion result = Quaternion()
        cdef Quaternion _p
        cdef Quaternion _q
        try:
            if isinstance(q, Quaternion) and isinstance(p, Quaternion):
                _p = p
                _q = q
                mult_quat_quat(&_q._components, &_p._components,
                        &result._components)
                return result
            if isinstance(q, Quaternion):
                _q = q
                mult_quat_scalar(&_q._components, p, &result._components)
                return result
            if isinstance(p, Quaternion):
                _p = p
                mult_quat_scalar(&_p._components, q, &result._components)
                return result
        except:
            pass

        return NotImplemented

    def __truediv__(Quaternion self, double scalar):
        return Quaternion(self.w/scalar, self.x/scalar, self.y/scalar, self.z/scalar)

    def __iadd__(Quaternion self, Quaternion q):
        self.w+=q.w
        self.x+=q.x
        self.y+=q.y
        self.z+=q.z
        return self

    def __isub__(Quaternion self, Quaternion q):
        self.w-=q.w
        self.x-=q.x
        self.y-=q.y
        self.z-=q.z
        return self

    def __imul__(Quaternion self,p):
        cdef Quaternion _p
        if isinstance(p,Quaternion):
            _p = p
            mult_quat_quat(&self._components, &_p._components, &self._components)
        else:
            self.w*=p
            self.x*=p
            self.y*=p
            self.z*=p
        return self

    property magnitude:
        def __get__(Quaternion self):
            return sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)

    property conjugate:
        def __get__(Quaternion self):
            return Quaternion(self.w,-self.x,-self.y,-self.z)

    property vector:
        def __get__(Quaternion self):
            cdef np.ndarray[np.float64_t, ndim=2] v = np.empty((3,1))
            cdef double* data = <double *> v.data
            data[0] = self._components.x
            data[1] = self._components.y
            data[2] = self._components.z
            return v

    def log(Quaternion self):
        """
        Natural logarithm of this quaternion.
        """
        cdef Quaternion result = Quaternion()
        quaternion_log(&self._components, &result._components)
        return result

    def exp(Quaternion self):
        """
        Exponential of this quaternion.
        """
        cdef Quaternion result = Quaternion()
        quaternion_exp(&self._components, &result._components)
        return result

    def __pow__(Quaternion q, double p,z):
        if p == -1:
            return q.conjugate/q.magnitude**2
        else:
            return (q.log()*p).exp()

    def copy(Quaternion self):
        """
        Obtain a copy of this quaternion.
        """
        return Quaternion(self.w,self.x,self.y,self.z,)

    __copy__ = copy

    def normalise(Quaternion self):
        """
        Normalise this quaternion to unit length.
        """
        scale = 1/self.magnitude
        self.w *= scale
        self.x *= scale
        self.y *= scale
        self.z *= scale
        return self

    def negate(Quaternion self):
        """
        Negate this quaternion.
        """
        self.w *= -1
        self.x *= -1
        self.y *= -1
        self.z *= -1

    def dot(Quaternion self, Quaternion other):
        """
        Compute the dot product of this quaternion with another.
        """
        return self.w*other.w + self.x*other.x + self.y*other.y + self.z*other.z

    def rotateVector(Quaternion self, np.ndarray[np.float64_t, ndim=2] v):
        """
        Rotate vectors using this quaternion.

        Equivalent to M{q*v*q.conjugate}.

        @param v: 3xN L{np.ndarray} of column vectors.

        @return: The vectors as they appear in their current co-ordinate frame
            after applying the rotation specified by this quaternion to them.
        """
        cdef np.ndarray[np.float64_t, ndim=2] result = np.empty((3,v.shape[1]))
        cdef quaternion_t pure
        cdef quaternion_t *p = &pure
        cdef quaternion_t *q = &self._components
        cdef quaternion_t *qconj = &(<Quaternion> self.conjugate)._components
        cdef int i

        for i in range(v.shape[1]):
            pure.w = 0
            pure.x = v[0, i]
            pure.y = v[1, i]
            pure.z = v[2, i]
            mult_quat_quat(q, p, p)
            mult_quat_quat(p, qconj, p)
            result[0, i] = pure.x
            result[1, i] = pure.y
            result[2, i] = pure.z

        return result

    def rotateFrame(Quaternion self, np.ndarray[np.float64_t, ndim=2] v):
        """
        Rotate co-ordinate frame of vectors using this quaternion.

        Equivalent to M{q.conjugate*v*q}.

        @param v: 3xN L{np.ndarray} of column vectors.

        @return: The vectors as they appear in the rotated co-ordinate frame
            obtained by applying the rotation specified by this quaternion to
            their current co-ordinate frame.
        """
        cdef np.ndarray[np.float64_t, ndim=2] result = np.empty((3,v.shape[1]))
        cdef quaternion_t pure
        cdef quaternion_t *p = &pure
        cdef quaternion_t *q = &self._components
        cdef quaternion_t *qconj = &(<Quaternion> self.conjugate)._components
        cdef int i

        for i in range(v.shape[1]):
            pure.w = 0
            pure.x = v[0, i]
            pure.y = v[1, i]
            pure.z = v[2, i]
            mult_quat_quat(qconj, p, p)
            mult_quat_quat(p, q, p)
            result[0, i] = pure.x
            result[1, i] = pure.y
            result[2, i] = pure.z

        return result

    def toMatrix(Quaternion self):
        """
        Obtain a 3x3 rotation matrix equivalent to the rotation specified by
        this quaternion.
        """
        m = np.matrix([[0,0,0],[0,0,0],[0,0,0]],dtype=float)
        m[0,0] = 1. - 2.*self.y**2 - 2.*self.z**2
        m[1,0] = 2. * (self.x*self.y + self.w*self.z)
        m[2,0] = 2. * (self.x*self.z - self.w*self.y)

        m[0,1] = 2. * (self.x*self.y - self.w*self.z)
        m[1,1] = 1. - 2.*self.x**2 - 2. *self.z**2
        m[2,1] = 2. * (self.y*self.z + self.w*self.x)

        m[0,2] = 2. * (self.x*self.z + self.w*self.y)
        m[1,2] = 2. * (self.y*self.z - self.w*self.x)
        m[2,2] = 1. - 2.*self.x**2 - 2.*self.y**2

        return m

    def toAxisAngle(Quaternion self):
        """
        Obtain the axis and angle of rotation specified by this quaternion.
        """
        angle = acos(self.w)
        axis = np.array((self.x,self.y,self.z))/sin(angle)
        return axis,2*angle

    def toEuler(Quaternion self, order='zyx',inDegrees=True):
        """
        Convert this quaternion to a corresponding Euler angle sequence.

        @param order: The order of the desired Euler angle sequence. Choices
            are 'zyx' for standard aerospace sequence (default),
            or 'zxy' for order used in BVH files.
        @param inDegrees: True to obtain angles in degrees (default), or False
            for radians.

        @return: Sequence of Euler angles.
        """
        m = self.toMatrix()
        return matrices.matrixToEuler(m,order,inDegrees)

    def set(Quaternion self, Quaternion other):
        """
        Set the compnents of this quaternion from another quaternion.
        """
        self.w = other.w
        self.x = other.x
        self.y = other.y
        self.z = other.z

    def setFromMatrix(Quaternion self, m):
        """
        Set this quaternion to be eqivalent to a given 3x3 rotation matrix.
        """
        t = m[0,0]+m[1,1]+m[2,2]
        if t > 0:
            w2 = sqrt(t+1)
            self.w = w2/2
            self.x = (m[2,1]-m[1,2])/(2*w2)
            self.y = (m[0,2]-m[2,0])/(2*w2)
            self.z = (m[1,0]-m[0,1])/(2*w2)
        else:
            t = m[0,0]-m[1,1]-m[2,2]
            if t > 0:
                x2 = sqrt(t+1)
                self.w = (m[2,1]-m[1,2])/(2*x2)
                self.x = x2/2
                self.y = (m[1,0]+m[0,1])/(2*x2)
                self.z = (m[0,2]+m[2,0])/(2*x2)
            else:
                t = m[1,1]-m[0,0]-m[2,2]
                if t > 0:
                    y2 = sqrt(t+1)
                    self.w = (m[0,2]-m[2,0])/(2*y2)
                    self.x = (m[1,0]+m[0,1])/(2*y2)
                    self.y = y2/2
                    self.z = (m[1,2]+m[2,1])/(2*y2)
                else:
                    z2 = sqrt(m[2,2]-m[0,0]-m[1,1]+1)
                    self.w = (m[1,0]-m[0,1])/(2*z2)
                    self.x = (m[0,2]+m[2,0])/(2*z2)
                    self.y = (m[1,2]+m[2,1])/(2*z2)
                    self.z = z2/2

    def setFromVectors(Quaternion self,x,y,z):
        """
        Set this quaternion so that it acts as a rotation taking the e1,e2,e3
        basis vectors to x,y,z.
        """
        self.setFromMatrix(np.hstack((x,y,z)).T)

    def setFromEuler(Quaternion self,angles,order='zyx',inDegrees=True):
        """
        Set this quaternion from an Euler angle sequence.

        @param angles: Sequence of 3 Euler angles.
        @param order: The order to apply the Euler angle sequence.
        @param inDegrees: True to indicate that angles are in degrees (default),
            or False for radians.
        """
        if inDegrees:
            angles = [np.radians(angle) for angle in angles]
        self.set( functools.reduce(operator.mul, [Quaternion(**dict((('w',cos(angle/2.0)),
            (axis.lower(),sin(angle/2.0))))) for angle, axis in zip(angles, order)]))



    fromEuler = staticmethod(QuaternionFromEuler)
    fromMatrix = staticmethod(QuaternionFromMatrix)
    fromVectors = staticmethod(QuaternionFromVectors)
    fromAxisAngle = staticmethod(QuaternionFromAxisAngle)
    nan = staticmethod(QuaternionNaN)

    def __reduce__(Quaternion self):
        return (Quaternion, (self.w, self.x, self.y, self.z))

def QuaternionArrayNaN(length):
    a = np.empty((length,4))
    a[:] = np.nan
    return QuaternionArray(a)

class QuaternionArray(object):
    """
    An array of quaternion values.

    Math operators are overridden to support quaternion math operations.
    """
    __slots__ = ['array','w','x','y','z']
    __hash__ = None
    def __init__(self,data,copy=False):
        """
        Construct quaternion array.

        A L{QuaternionArray} can be constructed from:

            - A sequence of L{Quaternion} objects.
            - A 4-element list or tuple of arrays giving w,x,y,z components.
            - An Nx4 array of quaternion component values.
            - Another L{QuaternionArray}.

        @param data: Data from which to construct the array.
        @param copy: If the source is an existing Nx4 array or
            L{QuaternionArray},  whether to copy the data. If not it will be
            referenced in place.
        """
        if isinstance(data, (tuple,list)):
            if isinstance(data[0],Quaternion):
                self.array = np.array([q.components for q in data],ndmin=2)
            elif len(data) == 4:
                self.array = np.array(data,ndmin=2)
            else:
                raise TypeError("List or tuple data must be a sequence of \
Quaternions or a sequence of w,x,y,z arrays")

        elif isinstance(data, np.ndarray):
            if isinstance(data[0], Quaternion):
                self.array = np.array([q.components for q in data],ndmin=2)
            else:
                assert data.ndim == 2 and data.shape[1] == 4, \
                    "Array data must be Nx4, got %r" % repr(data)
                self.array = np.array(data,ndmin=2,copy=copy)
        elif isinstance(data,QuaternionArray):
            self.array = np.array(data.array,copy=copy)
        elif isinstance(data,Quaternion):
            self.array = np.array([data.components], ndmin=2)
        else:
            raise TypeError("Cannot construct QuaternionArray from type %r"%type(data))

        shape = self.array.shape
        if not shape[1] == 4:
            self.array = self.array.T

        self.w = self.array[:,0]
        self.x = self.array[:,1]
        self.y = self.array[:,2]
        self.z = self.array[:,3]

    def __repr__(self):
        return "QuaternionArray(%r)"%self.array

    def __add__(self,other):
        if isinstance(other,(QuaternionArray,Quaternion)):
            return QuaternionArray((self.w+other.w,
                self.x+other.x, self.y+other.y, self.z+other.z))
        elif np.isscalar(other):
            return QuaternionArray(self.array+other)
        else:
            raise NotImplementedError

    def __sub__(self,other):
        if isinstance(other,(QuaternionArray,Quaternion)):
            return QuaternionArray((self.w-other.w,
                self.x-other.x, self.y-other.y, self.z-other.z))
        elif np.isscalar(other):
            return QuaternionArray(self.array-other)
        else:
            raise NotImplementedError

    def __neg__(self):
        return QuaternionArray(-self.array)

    def __mul__(self,other):
        if isinstance(other,(QuaternionArray,Quaternion)):
            return QuaternionArray((
                self.w*other.w - self.x*other.x - self.y*other.y - self.z*other.z,
                self.w*other.x + self.x*other.w + self.y*other.z - self.z*other.y,
                self.w*other.y - self.x*other.z + self.y*other.w + self.z*other.x,
                self.w*other.z + self.x*other.y - self.y*other.x + self.z*other.w))
        elif np.isscalar(other):
            return QuaternionArray(other*self.array)
        elif isinstance(other,np.ndarray):
            if other.ndim == 1 and other.shape[0] == self.array.shape[0]:
                return QuaternionArray(self.array*other.reshape(-1,1))
        else:
            raise NotImplementedError

    def __rmul__(self,other):
        if np.isscalar(other):
            return QuaternionArray(other*self.array)
        elif isinstance(other,Quaternion):
            return QuaternionArray((
                other.w*self.w - other.x*self.x - other.y*self.y - other.z*self.z,
                other.w*self.x + other.x*self.w + other.y*self.z - other.z*self.y,
                other.w*self.y - other.x*self.z + other.y*self.w + other.z*self.x,
                other.w*self.z + other.x*self.y - other.y*self.x + other.z*self.w))

        else:
            raise NotImplementedError

    def copy(self):
        return QuaternionArray(self.array.copy())

    __copy__ = copy

    def log(self):
        """
        Natural logarithm of the quaternions.
        """
        v = self.vector
        normv = vectors.norm(v)
        m = self.norm
        w = np.log(m)
        x,y,z = np.nan_to_num((v/normv)*np.arccos(self.w/m))
        return QuaternionArray((w,x,y,z))

    def exp(self):
        """
        Exponential of the quaternions.
        """
        v = self.vector
        normv = vectors.norm(v)
        w = np.cos(normv)
        x,y,z = np.nan_to_num(v/normv)*np.sin(normv)
        return QuaternionArray((w,x,y,z))*np.exp(self.w)

    def __pow__(self,p):
        return (self.log()*p).exp()

    def __eq__(self,other):
        if isinstance(other,QuaternionArray):
            return np.allclose(self.array,other.array)
        else:
            return False

    def __getitem__(self,key):
        result = self.array[key]
        if np.shape(result) == (4,):
            return Quaternion(*result)
        elif len(np.shape(result)) == 2 and np.shape(result)[1] == 4:
            return QuaternionArray(result)
        else:
            return result

    def __setitem__(self, key, value):
        if isinstance(value, Quaternion):
            self.array[key] = value.components
        elif isinstance(value, QuaternionArray):
            self.array[key] = value.array
        else:
            self.array[key] = value

    def __len__(self):
        return len(self.array)

    def dot(self,other):
        """
        Dot product with another quaternion array.
        """
        assert isinstance(other,QuaternionArray), 'Only QuaternionArrays \
are supported for dot product'
        return np.sum(self.array*other.array,axis=1)

    @property
    def magnitude(self):
        """
        Magnitudes of the quaternions.
        """
        return np.sqrt(self.dot(self))

    def rotateVector(self,v):
        """
        Rotate vectors by the rotations of these quaternions.

        Equivalent to M{q*v*q.conjugate} for each q and v.

        @param v: 3xN L{np.ndarray} of column vectors.

        @return: The vectors as they appear in their current co-ordinate frame
            after applying the rotations specified by these quaternions to
            them.
        """
        x,y,z = v.reshape((3,-1))
        r = np.empty((3,self.array.shape[0]))
        W = -self.x * x - self.y * y - self.z * z;
        X = self.w * x + self.y * z - self.z * y;
        Y = self.w * y - self.x * z + self.z * x;
        Z = self.w * z + self.x * y - self.y * x;

        r[2] = -W * self.z - X * self.y + Y * self.x + Z * self.w;
        r[1] = -W * self.y + X * self.z + Y * self.w - Z * self.x;
        r[0] = -W * self.x + X * self.w - Y * self.z + Z * self.y;
        return r

    def rotateFrame(self,v):
        """
        Rotate co-ordinate frames of vectors using these quaternions.

        Equivalent to M{q.conjugate*v*q} for each q and v.

        @param v: 3xN L{np.ndarray} of column vectors, where N equals the
            length of this quaternion array.

        @return: The vectors as they appear in the rotated co-ordinate frames
            obtained by applying the rotations given by these quaternions to
            their current co-ordinate frame.
        """
        x,y,z = v.reshape((3,-1))
        r = np.empty((3,self.array.shape[0]))
        W = self.x * x + self.y * y + self.z * z
        X = self.w * x - self.y * z + self.z * y
        Y = self.w * y + self.x * z - self.z * x
        Z = self.w * z - self.x * y + self.y * x

        r[2] = W * self.z + X * self.y - Y * self.x + Z * self.w
        r[1] = W * self.y - X * self.z + Y * self.w + Z * self.x
        r[0] = W * self.x + X * self.w + Y * self.z - Z * self.y
        return r

    @property
    def vector(self):
        """
        The (imaginary) column vector components of these quaternions.
        """
        return self.array[:,1:].T
    @property
    def conjugate(self):
        """
        The conjugates of these quaternions.
        """
        return QuaternionArray((self.w,-self.x,-self.y,-self.z))

    @property
    def norm(self):
        """
        The 2-norms of these quaternions.
        """
        return np.sqrt(np.sum(self.array**2,axis=1))

    def validity(self):
        """
        Obtain a boolean array indicating which elements are valid.
        """
        return ~np.any(np.isnan(self.array),axis=1)

    def unflipped(self):
        """
        Obtain a copy of this array with no sign flips between values.
        """
        new = QuaternionArray(self.array.copy())
        valid = self.validity()
        vals = self[valid]
        idx = np.flatnonzero(valid)
        for i in range(1,len(idx)):
            new[idx[i]] = -vals[i] if vals[i].dot(new[idx[i-1]]) < 0 else vals[i]
        new[~valid] = np.nan
        return new

    def smoothed(self, stddev=0.001):
        """
        Obtain a smoothed array of quaternions.

        Smoothing is perfomed using a cubic spline for each quaternion
        component. It is assumed that the time between quaternions is
        constant.

        @param stddev: Standard deviation of the noise expected in each
            component.
        """

        smoothedData = np.empty_like(self.array)
        ts = range(len(self.array))
        spl = PartialInputVectorSpline(ts, self.array.T, stddev=stddev)
        qa = QuaternionArray(spl(ts).T)
        return qa * (1.0/qa.magnitude)

    def toAxisAngle(self):
        """
        Obtain the axis and angle of rotations specified by these quaternions.
        """
        angle = np.arccos(self.w)
        axis = np.array((self.x,self.y,self.z))/np.sin(angle)
        return axis,2*angle

    nan = staticmethod(QuaternionArrayNaN)


cpdef QuaternionFactory(w,x,y,z):
    """
    Factory to create Quaternions or QuaternionArrays depending on the length
    of the arguments.
    """

    try:
        return Quaternion(w,x,y,z)
    except TypeError:
        if np.isscalar(w):
            w = np.zeros_like(x)
        return QuaternionArray(np.vstack((w,x,y,z)))
