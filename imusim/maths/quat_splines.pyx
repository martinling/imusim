# cython: profile=True
"""
Spline fitting of quaternion data.
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
from imusim.maths.quaternions import QuaternionArray, QuaternionFactory
from imusim.maths.splines import Spline, PartialInputSpline
import numpy as np
import math

from imusim.maths.quaternions cimport Quaternion, quaternion_t, mult_quat_quat
from imusim.maths.quaternions cimport mult_quat_scalar, quaternion_add, quaternion_exp
cimport numpy as np
cimport cython

DEF SPLINE_ORDER = 4

cdef class QuaternionBSpline:
    """
    Model of a quaternion function of time using B-spline fitting of keyframes.

    Quaternion interpolation is performed using the algorithm from "A general
    Construction Scheme for Unit Quaternion Curves with Simple High Order
    Derivatives" by Kim, M.-J., Kim M.-S., and Shin, S. Y.

    This interpolation method should be used in preference to the SQUAD
    algorithm as it is C^2 continuous.
    """
    cdef double tmin
    cdef double dt
    cdef list w
    cdef list q
    cdef double _validFrom
    cdef double _validTo

    def __cinit__(self, timestamps, quaternions, w=None):
        """
        Construct quaternion spline.

        @param timestamps: Sequence of monotonically increasing keyframe times.
        @param quaternions: L{QuaternionArray} of keyframe quaternions.
        """

        if len(timestamps) < 5:
            raise Spline.InsufficientPointsError

        self.tmin = timestamps[0]
        self.dt = timestamps[1] - timestamps[0]

        if w is None:
            self.q = [ q for q in quaternions ]
            # Extend quaternions to provide control points at the start and end
            self.q.extend([quaternions[-1]]*2)
            self.q.extend([quaternions[0]]*2)

            self.w = [(self.q[i-1]**-1 * self.q[i]).log() for i in
                    range(len(self.q))]
        else:
            self.q = quaternions
            self.w = w

        self._validFrom = timestamps[2]
        self._validTo = timestamps[-2]

    property validFrom:
        def __get__(self):
            return self._validFrom

    property validTo:
        def __get__(self):
            return self._validTo

    def __reduce__(self):
        return (QuaternionBSpline, ([self.tmin, self.tmin + self.dt], self.q, self.w))

    @cython.profile(False)
    cdef double B(self, int i,int k,double t):
        """
        B-Spline basis function B_{i,k}(t).
        """
        if k == 1:
            return 1 if i <= t < i+1 else 0

        return ((t-i)/(k-1)) * self.B(i,k-1,t) + \
                ((i+k-t)/(k-1)) * self.B(i+1,k-1,t)

    @cython.profile(False)
    cdef inline double Bprime(self, int i, int k, double t):
        """
        1st derivative of basis function B_{i,k}(t).
        """
        return (1/self.dt)*(self.B(i,k-1,t) - self.B(i+1,k-1,t))

    @cython.profile(False)
    cdef inline double tildeB(self, int i,int k,double t):
        """
        Cumulative basis function \\tilde{B}_{i,k}(t).
        """
        if t >= i+k-1:
            return 1
        if t <= i:
            return 0

        cdef double sum = 0.0
        cdef int j
        for j in range(i,i+k+1):
            sum += self.B(j,k,t)
        return sum

    def __call__(self, t):
        """
        Evaluate 0th, 1st and 2nd derivatives of the spline at times t.

        @return: L{QuaternionArray} of orientations, 3xN L{np.ndarray} of
            angular rate vectorsm 3xN L{np.ndarray} of angular acceleration
            vectors.
        """
        t -= self.tmin
        if np.isscalar(t):
            q,w,a = self.evaluate(t/self.dt)
        else:
            q = QuaternionArray(np.empty((len(t),4)))
            w = QuaternionArray(np.empty((len(t),4)))
            a = QuaternionArray(np.empty((len(t),4)))
            for i,T in enumerate(t/self.dt):
                q[i],w[i],a[i] = self.evaluate(T)

        return q, w.vector, a.vector

    cdef inline evaluate(self, double t):
        """
        Evaluate the derivatives of this QuaternionBSpline at time t.
        """
        t += 2
        cdef int l = int(math.floor(t))
        cdef int i
        cdef int j
        cdef double B
        cdef Quaternion q0 = self.q[l-SPLINE_ORDER+1]
        cdef Quaternion q = q0.copy()
        cdef quaternion_t *w
        cdef quaternion_t exp[SPLINE_ORDER]
        cdef quaternion_t dot_exp[SPLINE_ORDER]
        cdef quaternion_t ddot_exp[SPLINE_ORDER]
        cdef quaternion_t wexp
        cdef quaternion_t tmp

        for j in range(SPLINE_ORDER-1):
            i = l - SPLINE_ORDER + 2 + j
            w = &(<Quaternion>self.w[i])._components
            B = self.B(i, SPLINE_ORDER-1, t)

            #exp[j] = (w * self.tildeB(i, SPLINE_ORDER, t)).exp()
            mult_quat_scalar(w, self.tildeB(i, SPLINE_ORDER, t), &exp[j])
            quaternion_exp(&exp[j], &exp[j])

            #wexp = w * exp[j]
            mult_quat_quat(w, &exp[j], &wexp)

            #dot_exp[j] = wexp * B
            mult_quat_scalar(&wexp, B, &dot_exp[j])

            #ddot_exp[j] = wexp * self.Bprime(i, SPLINE_ORDER-1, t) +
            #    w * dot_exp[j] * B
            mult_quat_scalar(&wexp, self.Bprime(i, SPLINE_ORDER-1, t),
                    &ddot_exp[j])
            mult_quat_quat(w, &dot_exp[j], &tmp)
            mult_quat_scalar(&tmp, B, &tmp)
            quaternion_add(&ddot_exp[j], &tmp, &ddot_exp[j])

            # q *= exp[j]
            mult_quat_quat(&q._components, &exp[j], &q._components)

        #-----------------------------------------------------------------------
        # \dot{q}
        cdef quaternion_t dotq
        mult_quat_quat(&dot_exp[0], &exp[1], &dotq)
        mult_quat_quat(&dotq, &exp[2], &dotq)

        mult_quat_quat(&exp[0], &dot_exp[1], &tmp)
        mult_quat_quat(&tmp, &exp[2], &tmp)
        quaternion_add(&dotq, &tmp, &dotq)

        mult_quat_quat(&exp[0], &exp[1], &tmp)
        mult_quat_quat(&tmp, &dot_exp[2], &tmp)
        quaternion_add(&dotq, &tmp, &dotq)

        mult_quat_quat(&q0._components, &dotq, &dotq)
        mult_quat_scalar(&dotq, (1.0 / self.dt), &dotq)

        #----------------------------------------------------------------------
        # \ddot{q}
        cdef quaternion_t ddotq_1
        mult_quat_quat(&ddot_exp[0], &exp[1], &ddotq_1)
        mult_quat_quat(&ddotq_1, &exp[2], &ddotq_1)

        mult_quat_quat(&exp[0], &ddot_exp[1], &tmp)
        mult_quat_quat(&tmp, &exp[2], &tmp)
        quaternion_add(&ddotq_1, &tmp, &ddotq_1)

        mult_quat_quat(&exp[0], &exp[1], &tmp)
        mult_quat_quat(&tmp, &ddot_exp[2], &tmp)
        quaternion_add(&ddotq_1, &tmp, &ddotq_1)

        mult_quat_quat(&q0._components, &ddotq_1, &ddotq_1)
        mult_quat_scalar(&ddotq_1, (1.0 / self.dt), &ddotq_1)

        cdef quaternion_t ddotq_2
        mult_quat_quat(&dot_exp[0], &dot_exp[1], &ddotq_2)
        mult_quat_quat(&ddotq_2, &exp[2], &ddotq_2)

        mult_quat_quat(&dot_exp[0], &exp[1], &tmp)
        mult_quat_quat(&tmp, &dot_exp[2], &tmp)
        quaternion_add(&ddotq_2, &tmp, &ddotq_2)

        mult_quat_quat(&exp[0], &dot_exp[1], &tmp)
        mult_quat_quat(&tmp, &dot_exp[2], &tmp)
        quaternion_add(&ddotq_2, &tmp, &ddotq_2)

        mult_quat_quat(&q0._components, &ddotq_2, &ddotq_2)
        mult_quat_scalar(&ddotq_2, (2.0 / self.dt), &ddotq_2)

        quaternion_add(&ddotq_1, &ddotq_2, &ddotq_1)

        #-----------------------------------------------------------------------
        cdef quaternion_t qc2 = (<Quaternion> q.conjugate)._components
        mult_quat_scalar(&qc2, 2, &qc2)

        mult_quat_quat(&qc2, &dotq, &dotq)
        dq = Quaternion(dotq.w, dotq.x, dotq.y, dotq.z)

        mult_quat_quat(&qc2, &ddotq_1, &ddotq_1)
        ddq = Quaternion(ddotq_1.w, ddotq_1.x, ddotq_1.y, ddotq_1.z)

        return q, dq, ddq

class SQUADQuaternionSpline(Spline):
    """
    Model of a quaternion function of time using the SQUAD algorithm.

    The implementation is based on the derivation given in "Quaternion Algebra
    and Calculus" by David Eberly of Geometric Tools Inc.
    """

    def __init__(self, timestamps, quaternions):
        """
        Construct SQUAD interpolator.

        @param timestamps: Sequence of monotonically increasing keyframe times.
        @param quaternions: L{QuaternionArray} of keyframe quaternions.
        """
        assert len(timestamps) == len(quaternions)
        self._quaternions = [q for q in quaternions]
        self._timestamps = np.array(timestamps,copy=True)

        # Repeat start and end quaternions to allow interpolation in
        # first and last intervals. Result is that q[-1] = q[0] and
        # q[N] = q[N-1], where N is the number of interpolation points
        self._quaternions.append(self._quaternions[-1])
        self._quaternions.append(self._quaternions[0])

        N = len(quaternions)
        assert self._quaternions[0] == self._quaternions[-1]
        assert self._quaternions[N] == self._quaternions[N-1]

        self._controlPoints = self._calculateControlPoints()

        dts = np.diff(timestamps)
        self._dt = dts[0] if np.allclose(dts,dts[0]) else None

    def _indexfunction(self,t):
        if self._dt:
            return int(math.ceil(t/self._dt))
        else:
            return np.searchsorted(self._timestamps,t)

    def __call__(self,t,n):
        """
        Evaluate the n-th derivative of the quaternion function at times t.
        """
        if n == 0:
            return self.squad(*self._setupParams(t))
        elif n == 1:
            return self.squadOmega(*self._setupParams(t))
        else:
            raise NotImplementedError('Order %d derivatives not implemented')

    def _setupParams(self,ts):
        """
        Setup parameters for squad and squadPrime calls.
        """
        ts = np.array(ts,ndmin=1)
        shape = (ts.shape[0],4)
        q1 = np.empty(shape)
        q2 = np.empty(shape)
        a = np.empty(shape)
        b = np.empty(shape)
        u = np.empty_like(ts)
        dt = np.empty_like(ts)
        for i,t in enumerate(ts):
            futureIndex = self._indexfunction(t)
            pastIndex = futureIndex-1
            q1[i] = self._quaternions[futureIndex].components
            q2[i] = self._quaternions[pastIndex].components
            a[i] = self._controlPoints[pastIndex].components
            b[i] = self._controlPoints[futureIndex].components
            u[i] = ((t-self._timestamps[pastIndex]) /
                (self._timestamps[futureIndex] -self._timestamps[pastIndex]))
            dt = self._timestamps[futureIndex] - self._timestamps[pastIndex]

        return (QuaternionArray(q1),QuaternionArray(q2),QuaternionArray(a),
            QuaternionArray(b),u,dt)

    @staticmethod
    def slerp(q1,q2,u):
        """
        Perform spherical linear interpolation between quaternions q1 and q2 at
        normalised time u.

        @param q1: Starting quaternion.
        @param q2: Ending quaternion.
        @param u: Normalised time parameter 0 <= u <= 1.
        """
        assert np.all(u >= 0)
        assert np.all(u <= 1)

        cosTheta = q1.dot(q2)
        theta = np.arccos(cosTheta)
        sinTheta = np.sin(theta)
        return (q1*np.nan_to_num(np.sin((1-u)*theta)/sinTheta) +
            q2*np.nan_to_num(np.sin(u*theta)/sinTheta))

    @staticmethod
    def slerpPrime(q1,q2,u):
        """
        Calculate the derivative of the slerp between q1 and q2 at normalised
        time u.
        """
        q = q1**-1 * q2
        return q1*q**u * q.log()

    @staticmethod
    def slerpOmega(q1,q2,u):
        """
        Calculate the rotational rate vector, in the rotating frame, at
        normalised time u as slerping between q1 and q2.

        Uses the relationship M{qdot = 1/2 omega * q}.
        """
        qdot = SQUADQuaternionSpline.slerpPrime(q1,q2,u)
        q = SQUADQuaternionSpline.slerp(q1,q2,u)
        return (2 * qdot*q.conjugate).vector

    @staticmethod
    def squad(q1,q2,a,b,u,*args):
        """
        Perform spherical cubic interpolation.

        @param q1: Starting quaternion.
        @param q2: Ending quaternion
        @param a: Control point a.
        @param b: Control point b.
        @param u - normalised time 0 <= u <= 1
        @param args: Additional arguments are ignored.
        """
        slerp = SQUADQuaternionSpline.slerp
        return slerp(slerp(q1,q2,u),slerp(a,b,u),2*u*(1-u))

    @staticmethod
    def squadPrime(q1,q2,a,b,u):
        """
        Calculate the derivative of SQUAD at normalised time u.
        """
        slerp = SQUADQuaternionSpline.slerp
        slerpPrime = SQUADQuaternionSpline.slerpPrime
        U = slerp(q1,q2,u)
        V = slerp(a,b,u)
        W = U**-1 * V
        UPrime = slerpPrime(q1,q2,u)
        VPrime = slerpPrime(a,b,u)
        WPrime = U**-1 * VPrime - U**-2 * UPrime *V
        t = 2*u*(1-u)
        return U*(W**t*W.log()*(2-4*u) + W**(t-1)*WPrime*t) + UPrime*W**t

    @staticmethod
    def squadOmega(q1,q2,a,b,u,dt):
        """
        Calculate the rotational rate vector, in the rotating frame, at
        normalised time u, while performing SQUAD(q1,q2,a,b,u).
        """
        qdot = SQUADQuaternionSpline.squadPrime(q1,q2,a,b,u)
        q = SQUADQuaternionSpline.squad(q1,q2,a,b,u)
        return (2 * qdot*q.conjugate).vector/dt

    def _calculateControlPoints(self):
        """
        Calculate the appropriate control point for Squad interpolation given
        the index, n, into the keyframe array.
        """
        q = self._quaternions
        N = len(q)-1

        return [q[n] * np.exp(-(np.log(q[n]**-1 * q[n+1]) +
            np.log(q[n]**-1 * q[n-1]))/4) for n in range(N)]

class PartialInputQuaternionBSpline(PartialInputSpline):
    """
    Quaternion B-spline allowing for undefined regions in output domain.

    Call the resulting object returns a tuple of a QuaternionArray,
    rotational rate array and rotational acceleration array.
    """
    _splineClass = QuaternionBSpline

    def _validity(self, y):
        return y.validity()

    def _output(self, x, conditions, results, undefined):
        length = len(np.atleast_1d(x))
        q = QuaternionArray(np.empty((length,4)))
        w = np.empty((3,length))
        a = np.empty((3,length))
        for condition, result in zip(conditions, results):
            q[condition], w[:,condition], a[:,condition] = result
        q.array[undefined] = np.nan
        w[:,undefined] = np.nan
        a[:,undefined] = np.nan
        if np.isscalar(x):
            return q[0],w,a
        else:
            return q,w,a
