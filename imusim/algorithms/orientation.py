"""
Algorithms for tracking orientation using inertial/magnetic sensor data.
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

from abc import ABC, abstractmethod
from imusim.maths.quaternions import Quaternion
from imusim.maths import vectors
from imusim.maths.kalman import KalmanFilter
from imusim.algorithms import vector_observation
from imusim.utilities.time_series import TimeSeries
from imusim.environment.gravity import STANDARD_GRAVITY
from imusim.utilities.documentation import prepend_method_doc
import copy
import collections
import numpy as np
import math


class OrientationFilter(ABC):
    """
    Base class for orientation estimation filters.

    An orientation filter takes accelerometer, magnetometer and gyroscope
    readings from an IMU and attempts to track the orientation of the IMU
    as a L{Quaternion} from these measurements.

    @ivar rotation: L{TimeSeries} of quaternion orientation estimates.
    """

    def __init__(self, initialTime, initialRotation):
        """
        Initialise orientation filter.

        @param initialTime: Initial time.
        @param initialRotation: Initial rotation (L{Quaternion}).
        """
        assert isinstance(initialRotation, Quaternion)
        self.rotation = TimeSeries()

    def __call__(self, accel, mag, gyro, t):
        """
        Run filter update.

        @param accel: Accelerometer sample (3x1 L{np.ndarray})
        @param mag: Magnetometer sample (3x1 L{np.ndarray})
        @param gyro: Gyroscope sample (3x1 L{np.ndarray})
        @param t: Time at which samples were taken.
        """
        dt = t - self.rotation.latestTime
        self._update(accel, mag, gyro, dt, t)

    @abstractmethod
    def _update(self, accel, mag, gyro, dt, t):
        """
        Internal filter update routine.

        Should add its results to the self.rotation L{TimeSeries}.
        """
        pass

class GyroIntegrator(OrientationFilter):
    """
    An estimator that obtains its orientation using only gyroscope integration.
    """

    def __init__(self, initialTime, initialRotation, **kwargs):
        OrientationFilter.__init__(self, initialTime, initialRotation)
        self.rotation.add(initialTime, initialRotation)
        self.qHat = initialRotation.copy()

    def _update(self, accel, mag, gyro, dt, t):
        dotq = 0.5 * self.qHat * Quaternion(0, *gyro)
        self.qHat += dotq * dt
        self.qHat.normalise()
        self.rotation.add(t, self.qHat)

class BachmannCF(OrientationFilter):
    """
    Implementation of the complementary filter algorithm proposed by Bachmann.

    See E. R. Bachmann. "Inertial and Magnetic Tracking of Limb Segment
    Orientation for Inserting Humans into Synthetic Environments". PhD thesis,
    Naval Postgraduate School, California, 2000.

    @ivar vectorObservation: L{TimeSeries} of vector observation results.
    """

    @prepend_method_doc(OrientationFilter)
    def __init__(self, initialTime, initialRotation, k, **kwargs):
        """
        @param k: Filter blending co-efficient (float). At each timestep, a
            correction will be applied of M{1/k} times the difference between
            the gyro integration and vector observation results.
        """
        OrientationFilter.__init__(self, initialTime, initialRotation)
        self.rotation.add(initialTime, initialRotation)
        self._k = float(k)
        self.qHat = initialRotation.copy()
        self._vectorObservation = vector_observation.FQA()
        self.vectorObservation = TimeSeries()

    def _update(self, accel, mag, gyro, dt, t):
        dotq = 0.5 * self.qHat * Quaternion(0, *gyro)
        self.qHat += dotq * dt
        qMeas = self._vectorObservation(-accel, mag)
        if self.qHat.dot(qMeas) < 0:
            qMeas.negate()
        qError = qMeas - self.qHat
        self.qHat += (1/self._k) * dt * qError
        self.qHat.normalise()
        self.vectorObservation.add(t, qMeas)
        self.rotation.add(t, self.qHat)

class OrientCF(OrientationFilter):
    """
    Implementation of the complementary filter used on the Orient IMU.

    See A. Young, M. Ling, and D. K. Arvind. "Orient-2: A Realtime Wireless
    Posture Tracking System using Local Orientation Estimation". in Proc.
    4th Workshop on Embedded Network Sensors, pp 53-57. ACM, 2007.

    @ivar vectorObservation: L{TimeSeries} of vector observation results.
    """

    @prepend_method_doc(BachmannCF)
    def __init__(self, initialTime, initialRotation, k, aT, **kwargs):
        """
        @param aT: acceleration threshold (float). The correction step will
            only be performed if the condition M{abs(norm(accel) - 1) <= aT}
            is met.
        """
        OrientationFilter.__init__(self, initialTime, initialRotation)
        self.rotation.add(initialTime, initialRotation)
        self._k = float(k)
        self._aT = float(aT)
        self.qHat = initialRotation.copy()
        self._vectorObservation = vector_observation.GramSchmidt()
        self.vectorObservation = TimeSeries()

    def _update(self, accel, mag, gyro, dt, t):
        dotq = 0.5 * self.qHat * Quaternion(0,*gyro)
        self.qHat += dotq * dt

        if abs(vectors.norm(accel) - 1) < self._aT:
            qMeas = self._vectorObservation(-accel, mag)
            if self.qHat.dot(qMeas) < 0:
                qMeas.negate()
            qError = qMeas - self.qHat
            self.qHat += (1/self._k) * dt * qError
        else:
            qMeas = Quaternion.nan()
        self.vectorObservation.add(t, qMeas)

        self.qHat.normalise()
        self.rotation.add(t, self.qHat)

class YunEKF(OrientationFilter):
    """
    Implementation of the extented Kalman filter proposed by Yun et al.

    See X. Yun and E. R. Bachmann. "Design, Implementation, and Experiemental
    Results of a Quaternion-Based Kalman Filter for Human Body Motion
    Tracking". IEEE Transactions on Robotics, vol 22, no 6, pp 1216-1227, 2006.

    @ivar vectorObservation: L{TimeSeries} of vector observation results.
    @ivar rotationalVelocity: L{TimeSeries} of angular rate estimates.
    """

    @prepend_method_doc(OrientationFilter)
    def __init__(self, initialTime, initialRotation, initialRotationalVelocity,
            initialCovariance, measurementCovariance, D, tau, **kwargs):
        """
        @param initialRotationalVelocity: Initial angular rate estimate
            (3x1 L{np.ndarray}).
        @param initialCovariance: Initial state covariance matrix (7x7).
        @param measurementCovariance: Measurement noise covariance (7x7).
        @param D: Variance of process noise model in rad^2/s^2 (float).
        @param tau: Time constant of process noise model in seconds (float).
        """
        OrientationFilter.__init__(self, initialTime, initialRotation)
        self.rotation.add(initialTime, initialRotation, initialCovariance[3:,3:])

        self.tau = float(tau)
        self.D = float(D)

        self._vectorObservation = vector_observation.FQA()

        # Process noise covariance matrix, updated on each call.
        self.Q = np.asmatrix(np.zeros((7,7)))

        # Measurement noise covariance matrix
        self.R = measurementCovariance

        # Prior covariance estimate and state vector initialisation.
        self.P_minus = initialCovariance
        self.xHat_minus = np.vstack((initialRotationalVelocity,
            np.vstack(initialRotation.components)))

        self.vectorObservation = TimeSeries()
        self.rotationalVelocity = TimeSeries()

    def _update(self, accel, mag, gyro, dt, t):
        # Note: Measurement matrix H is equivalent to the identity matrix
        # so is left out to simplify calculations

        # Calculate Kalman gain
        K = self.P_minus * np.linalg.inv(self.P_minus + self.R)

        # Construct measurment vector z
        self.vo = self._vectorObservation(-accel, mag)
        if self.vo.dot(Quaternion(*self.xHat_minus[3:])) < 0:
            self.vo.negate()
        z = np.vstack((gyro, np.vstack(self.vo.components)))

        # Calculate state update
        xHat = self.xHat_minus + K * (z - self.xHat_minus)
        P = (np.eye(7) - K) * self.P_minus

        # Normalise quaternion component of state vector
        qHat = Quaternion(*xHat[3:]).normalise()
        xHat[3:] = np.vstack(qHat.components)

        # Calculate prediction
        self.xHat_minus = xHat + self.dotx(xHat) * dt
        Phi = self.Phi(xHat, dt)
        self.Q[[0,1,2],[0,1,2]] = \
            (self.D / (2 * self.tau)) * (1 - math.e**((-2 * dt) / self.tau))
        self.P_minus = Phi * P * Phi.T + self.Q

        self.rotation.add(t, qHat, P[3:,3:])
        self.vectorObservation.add(t, self.vo)
        self.rotationalVelocity.add(t, xHat[:3], P[:3,:3])

    def dotx(self, xHat):
        """
        Calculate state derivative vector.
        """
        dotx = np.empty((7, 1))
        dotx[:3] = -xHat[:3] / self.tau
        dotx[3] = -0.5 * (xHat[0]*xHat[4] + xHat[1]*xHat[5] + xHat[2]*xHat[6])
        dotx[4] = 0.5 * (xHat[0]*xHat[3] - xHat[1]*xHat[6] +xHat[2]*xHat[5])
        dotx[5] = 0.5 * (xHat[0]*xHat[6] + xHat[1]*xHat[3] - xHat[2]*xHat[4])
        dotx[6] = 0.5 * (-xHat[0]*xHat[5] +xHat[1]*xHat[4] + xHat[2]*xHat[3])
        return dotx

    def Phi(self, xHat, dt):
        """
        Calculate discrete state transition matrix Phi.
        """
        # Note: Yun paper matrix element indices are one-based so subtract 1
        # from all for zero-based indices.
        Phi = np.eye(7)
        Phi[0,0] = math.e**(-dt / self.tau)
        Phi[1,1] = math.e**(-dt / self.tau)
        Phi[2,2] = math.e**(-dt / self.tau)

        Phi[3,0] = -xHat[4] * dt / 2
        Phi[3,1] = -xHat[5] * dt / 2
        Phi[3,2] = -xHat[6] * dt / 2
        Phi[3,4] = -xHat[0] * dt / 2
        Phi[3,5] = -xHat[1] * dt / 2
        Phi[3,6] = -xHat[2] * dt / 2

        Phi[4,0] = xHat[3] * dt / 2
        Phi[4,1] = -xHat[6] * dt / 2
        Phi[4,2] = xHat[5] * dt / 2
        Phi[4,2] = xHat[0] * dt / 2
        Phi[4,5] = xHat[2] * dt / 2
        Phi[4,6] = -xHat[1] * dt / 2

        Phi[5,0] = xHat[6] * dt / 2
        Phi[5,1] = xHat[3] * dt / 2
        Phi[5,2] = -xHat[4] * dt / 2
        Phi[5,3] = xHat[1] * dt / 2
        Phi[5,4] = -xHat[2] * dt / 2
        Phi[5,6] = xHat[0] * dt / 2

        Phi[6,0] = -xHat[5] * dt / 2
        Phi[6,1] = xHat[4] * dt / 2
        Phi[6,2] = xHat[3] * dt / 2
        Phi[6,3] = xHat[2] * dt / 2
        Phi[6,4] = xHat[1] * dt / 2
        Phi[6,5] = -xHat[0] * dt / 2

        return Phi

class DistLinAccelCF(OrientationFilter):
    """
    Implementation of the multi-IMU orientation algorithm by Young & Ling.

    See A. D. Young, M. J. Ling and D. K. Arvind. "Distributed Estimation of
    Linear Acceleration for Improved Accuracy in Wireless Inertial Motion
    Capture", in Proc. 9th ACM/IEEE International Conference on Information
    Processing in Sensor Networks, ACM, 2010, pp. 256-267.

    @ivar qMeas: Result of latest vector observation (L{Quaternion}).
    @ivar vectorObservation: L{TimeSeries} of vector observation results.
    @ivar jointAccel: Latest estimate of joint acceleration
        (3x1 L{np.ndarray}).
    @ivar jointAcceleration: L{TimeSeries} of acceleration estimates.
    """
    GRAVITY_VECTOR = vectors.vector(0, 0, STANDARD_GRAVITY)

    @prepend_method_doc(BachmannCF)
    def __init__(self, initialTime, initialRotation, k, joint, offset,
            **kwargs):
        """
        @param joint: The L{Joint} the IMU is attached to.
        @param offset: The offset of the IMU in the joint's co-ordinate frame
            (3x1 L{np.ndarray})
        """
        OrientationFilter.__init__(self, initialTime, initialRotation)
        self.rotation.add(initialTime, initialRotation)
        self.qHat = initialRotation.copy()
        self.k = k
        self.joint = joint
        self.offset = offset

        self._vectorObservation = vector_observation.GramSchmidt()
        self.qMeas = Quaternion.nan()
        self.jointAccel = vectors.nan()

        self.gyroFIFO = collections.deque([], maxlen=3)
        self.accelFIFO = collections.deque([], maxlen=2)
        self.magFIFO = collections.deque([], maxlen=2)

        self.isRoot = not joint.hasParent
        self.children = joint.childJoints

        self.vectorObservation = TimeSeries()
        self.jointAcceleration = TimeSeries()

    def handleLinearAcceleration(self, jointAcceleration, dt):
        """
        Perform drift correction based on incoming joint acceleration estimate.

        @param jointAcceleration: Acceleration estimate (3x1 L{np.ndarray}).
        """
        self.jointAccel = jointAcceleration

        if len(self.gyroFIFO) < 3:
            return

        # Estimate linear acceleration
        o = self.offset
        w = self.gyroFIFO
        q = self.qHat_minus_1
        a = (w[0] - w[2]) / (2 * dt)
        lt = vectors.cross(a, o)
        lr = vectors.dot(w[1], o) * w[1] - o * vectors.norm(w[1])**2
        l = (q.rotateFrame(self.jointAccel) + lt + lr)

        # Apply drift correction
        self.qMeas = self._vectorObservation(-(self.accelFIFO[1]- l),
                self.magFIFO[1])
        if q.dot(self.qMeas) < 0:
            self.qMeas.negate()
        q = q + 1.0/self.k * dt * (self.qMeas - q)
        dotq = 0.5 * self.qHat * Quaternion(0,*w[0])
        self.qHat = q + dotq * dt
        self.qHat.normalise()

    def childAcceleration(self, o, dt):
        """
        Compute acceleration for child joint.

        @param o: Offset of child joint (3x1 L{np.ndarray}).
        """
        w = self.gyroFIFO
        if len(w) < 3:
            return None
        q = self.qHat_minus_1
        a = (w[0] - w[2]) / (2 * dt)
        lt = vectors.cross(a, o)
        lr = vectors.dot(w[1], o) * w[1] - o * vectors.norm(w[1])**2
        l = self.jointAccel + q.rotateVector(lt + lr)
        return l

    def _update(self, accel, mag, gyro, dt, t):
        # Store measurements in queues for processing
        self.gyroFIFO.appendleft(gyro)
        self.accelFIFO.appendleft(accel)
        self.magFIFO.appendleft(mag)

        # save current qHat for drift correction
        self.qHat_minus_1 = copy.copy(self.qHat)

        # Inertial update
        dotq = 0.5 * self.qHat * Quaternion(0,*gyro)
        self.qHat += dotq * dt
        self.qHat.normalise()

        if self.isRoot and len(self.accelFIFO) >= 2:
            # Root can't receive acceleration data so must estimate it.
            measuredAccel = self.qHat_minus_1.rotateVector(self.accelFIFO[1])
            # Assume root joint acceleration was zero.
            self.jointAccel = np.zeros((3,1))
            # Calculate tangential acceleration for root joint IMU offset.
            linAccel = self.childAcceleration(self.offset, dt)
            if linAccel is not None:
                # Calculate expected measurement.
                expectedAccel = linAccel - self.GRAVITY_VECTOR
                # Set root joint acceleration to the residual.
                self.jointAccel = measuredAccel - expectedAccel
                self.handleLinearAcceleration(self.jointAccel, dt)

        self.rotation.add(t, self.qHat)
        self.vectorObservation.add(t, self.qMeas)
        self.jointAcceleration.add(t, self.jointAccel)
