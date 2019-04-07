"""
Algorithms for estimating orientation from reference vector observations.
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
from abc import ABC, abstractmethod
from imusim.maths.quaternions import Quaternion
from scipy.optimize import newton
import imusim.maths.vectors as vectors
from imusim.utilities.documentation import prepend_method_doc
import numpy as np
import math


class VectorObservation(ABC):
    """
    Base class for all vector observation methods.

    Vector observation methods calculate the L{Quaternion} that rotates a
    set of reference vectors to match the set of observed vectors.
    Vectors are assumed to be 3x1 L{np.ndarray}s
    """
    def __call__(self, *measurements):
        """
        Estimate the orientation Quaternion from the observed vectors.

        @param measurements: The observed vectors (3x1 L{np.ndarray})
        @return: The estimated orientation L{Quaternion}
        """
        if any(map(lambda a: np.any(np.isnan(a)), measurements)):
            return Quaternion(np.nan,np.nan,np.nan,np.nan)
        else:
            return self._process(*measurements)

    @abstractmethod
    def _process(self, *measurements):
        """
        Internal method to be implemented by subclasses

        Implementations should return the L{Quaternion} that rotates the
        reference vectors to match the passed vector observations

        @param measurements: The observed vectors (3x1 L{np.ndarray})
        @return: The estimated orientation L{Quaternion}
        """
        pass


class TRIAD(VectorObservation):
    """
    Implementation of the TRIAD vector observation algorithm.

    The algorithm is designed to estimation the orientation based on
    observations of the Earth gravitational and magnetic fields based on the
    assumption of a (X,Y,Z) -> (North,East,Down) co-ordinate frame

    See M.D. Shuster, and S.D. Oh "Three-Axis Attitude Determination from
    Vector Observations" in Journal of Guidance and Control v.4 pp 70-77.
    1981
    """

    def _process(self, g, m):
        """
        Estimate the orientation Quaternion from the observed vectors.

        @param g: the observed gravitational field vector (3x1 L{np.ndarray})
        @param m: the observed magnetic field vector (3x1 L{np.ndarray})
        @return: The estimated orientation L{Quaternion}
        """
        z = g / vectors.norm(g)
        y = vectors.cross(z, m)
        y /= vectors.norm(y)
        x = vectors.cross(y, z)

        return Quaternion.fromVectors(x, y, z)


class GramSchmidt(VectorObservation):
    """
    Vector observation based on Gram-Schmidt Ortho-Normalisation.

    The algorithm is designed to estimation the orientation based on
    observations of the Earth gravitational and magnetic fields based on the
    assumption of a (X,Y,Z) -> (North,East,Down) co-ordinate frame
    """

    @prepend_method_doc(TRIAD)
    def _process(self, g, m):
        z = g / vectors.norm(g)
        x = m - (z * vectors.dot(m, z))
        x /= vectors.norm(x)
        y = vectors.cross(z, x)

        return Quaternion.fromVectors(x, y, z)


class FQA(VectorObservation):
    """
    Implementation of the Factored Quaternion Algorithm by Yun et al.

    The algorithm is designed to estimation the orientation based on
    observations of the Earth gravitational and magnetic fields based on the
    assumption of a (X,Y,Z) -> (North,East,Down) co-ordinate frame

    See X. Yun, E. Bachmann, and R. McGhee "A Simplified Quaternion-Based
    Algorithm for Orientation Estimation from Earth Gravity and Magnetic
    Field Measurements" In IEEE Trans. on Instrumentation and Measurement
    v. 57 pp 638-650 2008
    """
    def __init__(self):
        """
        Initialise vector observation.
        """
        alpha = np.radians(20)
        self.qAlpha = Quaternion(math.cos(alpha/2), 0, math.sin(alpha/2), 0)
        self.thetaThreshold = 0.1

    def cosHalfAngle(self, cosAngle):
        """
        Compute the value of cos(angle/2) from cos(angle)

        @param cosAngle: cosine of angle (float)
        """

        if -1.0001 < cosAngle <= -1:
            cosAngle = -1
        return math.sqrt((1+cosAngle)/2)

    def sinHalfAngle(self, cosAngle, sinAngle):
        """
        Compute the value of sin(angle/2) from cos(angle) and sin(angle)

        @param cosAngle: cosine of angle (float)
        @param sinAngle: sine of angle (float)
        """
        if 1 <= cosAngle < 1.0001:
            cosAngle = 1
        sign =  1 if sinAngle >= 0 else -1
        return sign * math.sqrt((1-cosAngle)/2)

    @prepend_method_doc(TRIAD)
    def _process(self, g, m):
        # Algorithm assumes that gravity vector is -z in default position
        # This is contrary to the other algorithm so we invert
        g = -g
        g /= vectors.norm(g)
        # Pitch
        sinTheta = g[0]
        cosTheta = math.sqrt(1 - sinTheta**2)

        flag = cosTheta <= self.thetaThreshold
        if flag:
            # If flag is set then theta is close to zero. Rotate coord frame 20
            # degrees to reduce errors
            g = self.qAlpha.rotateFrame(g)
            m =  self.qAlpha.rotateFrame(m)

            sinTheta = g[0]
            cosTheta = math.sqrt(1 - sinTheta**2)

        cosHalfTheta = self.cosHalfAngle(cosTheta)
        sinHalfTheta = self.sinHalfAngle(cosTheta, sinTheta)
        qp = Quaternion(cosHalfTheta, 0, sinHalfTheta, 0)

        # Roll
        cosPhi =  -g[2] / cosTheta
        sinPhi =  -g[1] / cosTheta
        cosHalfPhi = self.cosHalfAngle(cosPhi)
        sinHalfPhi = self.sinHalfAngle(cosPhi, sinPhi)
        qr = Quaternion(cosHalfPhi, sinHalfPhi, 0, 0)
        # Yaw
        m = qr.rotateVector(m)
        m = qp.rotateVector(m)
        Mx = m[0]
        My = m[1]
        scale = 1 / math.sqrt(Mx**2 + My**2)
        Mx *= scale
        My *= scale

        # From Yun paper
        #cosPsi = Mx * self.Nx + My * self.Ny
        #sinPsi = Mx * self.Ny - My * self.Nx
        # but we know that Ny is 0 so this simplifies to
        cosPsi = Mx
        sinPsi = -My
        cosHalfPsi = self.cosHalfAngle(cosPsi)
        sinHalfPsi = self.sinHalfAngle(cosPsi, sinPsi)
        qy = Quaternion(cosHalfPsi, 0, 0, sinHalfPsi)
        combined = qy * qp * qr

        if flag:
            return combined * self.qAlpha.conjugate
        else:
            return combined


class LeastSquaresOptimalVectorObservation(VectorObservation):
    """
    Base class of least squares optimal vector observation algorithms

    Implementations provide optimal solutions to Wahba's loss function. See
    G. Wahba "A Least Squares Estimate of Satellite Attitude" SIAM review v.7
    pp 384-386 1965

    @ivar refs: L{np.ndarray} of reference vectors in world frame
    @ivar weights: L{np.ndarray} of weights to apply to reference vectors
    """

    def __init__(self, refs=None, weights=None, inclinationAngle=None):
        """
        Initialise vector observation algorithm.

        The set of reference vectors for the vector observation process can
        be supplied directly, or an magnetic field inclination angle can be
        supplied to indicate that the reference vectors should be the Earth
        gravitational and magnetic fields.

        An optional weights list can be supplied to support weighting of
        reference vectors when calculating the Wahba loss function.

        @param refs: list of reference vectors (3x1 L{np.ndarray})
        @param weights: list of reference vector wieghts (float)
        @param inclinationAngle: inclination angle of magnetic field (float)
        """
        assert refs is not None or inclinationAngle is not None, \
                "Must specify reference vectors or inclination angle"
        if refs is None:
            refs = [
                np.array([[0,0,1.0]]).T,
                np.array([[math.cos(math.radians(inclinationAngle)),
                    0, math.sin(math.radians(inclinationAngle))]]).T]
        if weights is None:
            weights = [1.0 for ref in refs]

        self.refs = np.asarray(refs)
        self.weights = np.asarray(weights) / sum(weights)


class DavenportQ(LeastSquaresOptimalVectorObservation):
    """
    Implementation of the Davenport-q algorithm.

    See F. Markley, and D. Mortari "Quaternion Attitude Estimation using
    Vector Observations" Journal of the Astronomical Sciences v.48 pp 359-380
    2000
    """
    def _process(self, *meas):
        """
        Estimate the orientation Quaternion from the observed vectors.

        @param meas: The observed vectors (3x1 L{np.ndarray})
        @return: The least squares optimal L{Quaternion}
        """
        B = sum([a * np.dot(b, r.T)
            for a,b,r in zip(self.weights, meas, self.refs)])
        S = B + B.T
        z = sum([a * vectors.cross(b, r)
            for a,b,r in zip(self.weights, meas, self.refs)])

        sigma = np.trace(B)
        K = np.empty((4,4))
        K[0:3,0:3] = S-sigma*np.identity(3)
        K[0:3,3] = z[:,0]
        K[3,0:3] = z[:,0]
        K[3,3] = sigma

        eigenValues,eigenVectors = np.linalg.eig(K)
        q = np.asarray(eigenVectors[:,np.argmax(eigenValues)])

        return Quaternion(q[3],q[0],q[1],q[2]).normalise()


class QUEST(LeastSquaresOptimalVectorObservation):
    """
    Implementation of the QUEST algorithm.

    See M.D. Shuster, and S.D. Oh "Three-Axis Attitude Determination from
    Vector Observations" in Journal of Guidance and Control v.4 pp 70-77.
    1981
    """
    def _apply(self, refs, meas):

        B = sum([a * np.dot(m, r.T) for a,m,r in zip(self.weights, meas, refs)])
        S = B + B.T
        Z = sum([a * vectors.cross(m, r) for a,m,r in zip(self.weights, meas, refs)])

        sigma = np.trace(B)
        kappa = ((S[1,1]*S[2,2] - S[1,2]*S[2,1]) +
                 (S[0,0]*S[2,2] - S[0,2]*S[2,0]) +
                 (S[0,0]*S[1,1] - S[0,1]*S[1,0]))
        # where did the above come from and why doesn't this work instead?:
        # kappa = np.trace(np.matrix(S).H)
        delta = np.linalg.det(S)

        if len(refs) == 2:
            # Closed form solution for lambdaMax, see Shuster 1979 eqs 72-3.
            cosTerm = vectors.dot(*refs)[0]*vectors.dot(*meas)[0] + \
                vectors.norm(vectors.cross(*refs))*vectors.norm(vectors.cross(*meas))
            lambdaMax = math.sqrt(np.sum(self.weights**2) + 2*np.product(self.weights)*cosTerm)
        else:
            # Iterate for lambdaMax, using substitutions from Shuster JAS1290 paper
            # for better numerical accuracy.
            a = sigma**2 - kappa
            b = sigma**2 + np.dot(Z.T, Z)[0,0]
            c = 8 * np.linalg.det(B)
            d = np.dot(Z.T, np.dot(np.dot(S,S), Z))[0,0]
            f = lambda l: (l**2 - a)*(l**2 - b) - c*l + c*sigma - d
            fprime = lambda l: 4*l**3 -2*(a+b)*l + c
            lambdaMax = newton(f, 1.0, fprime)

        alpha = lambdaMax**2 - sigma**2 + kappa
        beta = lambdaMax - sigma
        gamma = (lambdaMax + sigma)*alpha - delta
        X = np.dot(alpha*np.identity(3) + beta*S + np.dot(S,S), Z)

        return [gamma] + list(X[0:3,0])

    @prepend_method_doc(DavenportQ)
    def _process(self, *meas):

        meas = [m / vectors.norm(m) for m in meas]

        EPS = 1e-8

        # Rotate the frame as necessary to avoid singularities.
        # See Shuster 1979 eqs 74-76.

        # Note differences in quaternion permutations from Shuster's paper.
        # These ones work, why don't the ones given in the paper?
        p = self._apply(self.refs, meas)
        if abs(p[0]) < EPS:
            refs = [r*np.array([[1,-1,-1]]).T for r in self.refs]
            p = self._apply(refs, meas)
            if abs(p[0]) < EPS:
                refs = [r*np.array([[-1,1,-1]]).T for r in self.refs]
                p = self._apply(refs, meas)
                if abs(p[0]) < EPS:
                    refs = [r*np.array([[-1,-1,1]]).T for r in self.refs]
                    p = self._apply(refs, meas)
                    q = Quaternion(p[3], p[2], -p[1], -p[0])
                else:
                    q = Quaternion(p[2], -p[3], -p[0], p[1])
            else:
                q = Quaternion(-p[1], p[0], -p[3], p[2])
        else:
            q = Quaternion(p[0], p[1], p[2], p[3])

        q.normalise()

        return q
