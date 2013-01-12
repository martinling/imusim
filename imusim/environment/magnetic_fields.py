"""
Magnetic field models.
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
import imusim.maths.vectors as vectors
from imusim.maths.vector_fields import VectorField, ConstantVectorField
from imusim.maths.transforms import cartesianToPolar, polarToCartesian
import math
import numpy as np

class EarthMagneticField(ConstantVectorField):
    """
    Simple model of the Earth's magnetic field.

    The Earth field is modelled as a constant vector field configured by
    magnitude, inclination from the horizontal and declination from true
    North.
    """

    def __init__(self, magnitude=50e-6, inclination=70, declination=0):
        """
        Construct field model.

        Default parameters result in a field with a magnitude of 50
        microTesla, a 70 degree downwards inclination and 0 degrees
        declination. These defaults are based on the nominal values for
        Edinburgh, UK, with declination set to 0 degrees so that magnetic
        and true North are equivalent.

        @param magnitude: Magnitude of the field, in Teslas.
        @param inclination: Inclination of the field to the horizontal,
            in degrees.
        @param declination: Delcination of the field from true North,
            in degrees.
        """
        theta = math.radians(declination)
        ctheta = math.cos(theta)
        stheta = math.sin(theta)
        phi = math.radians(inclination)
        cphi = math.cos(phi)
        sphi = math.sin(phi)
        value = magnitude * vectors.vector(cphi*ctheta, cphi*stheta, sphi)
        ConstantVectorField.__init__(self, value)

class DistortedMagneticField(VectorField):
    """
    Magnetic field with distortions superimposed.
    """

    def __init__(self, baseField, distortions=[]):
        """
        Construct distorted field model.

        @param baseField: The base field before distortions.
        @param distortions: List of distortion fields to superimpose.
        """
        self._baseField = baseField
        self._distortions = distortions

    def addDistortion(self, distortion):
        self._distortions.append(distortion)

    def __call__(self, position, t):
        return sum([d(position, t) for d in self._distortions],
                self._baseField(position, t))

    @property
    def nominalValue(self):
        return self._baseField.nominalValue

    def headingVariation(self, position, t):
        """
        Calculate heading variation due to magnetic distortion.

        @param position: 3xN L{np.ndarray} of position(s).

        @return: Heading variation(s) in degrees at given positions.
        """
        Bx,By,Bz = self(position, t)

        hField = np.vstack((Bx.ravel(),By.ravel()))
        hField /= vectors.norm(hField)
        ref = np.empty_like(hField)
        ref[:] = self.nominalValue[:2].copy()
        ref /= vectors.norm(ref)
        variation = np.arccos(vectors.dot(hField,ref))
        return np.degrees(variation)

class SolenoidMagneticField(object):
    """
    Model of magnetic field around an ideal solenoid.

    Follows equations from N. Derby and S. Olbert, "Cylindrical mangets and
    ideal solenoids", American Journal of Physics, vol 78, no. 3, pp. 229-235,
    March 2010.
    """
    def __init__(self,n,I,a,b,transform):
        """
        Construct solenoid field model.

        @param n: Number of coils per unit length.
        @param I: Coil current.
        @param a: Radius of the coil
        @param b: Length of the solenoid.
        @param transform: L{AffineTransform} specifying the orientation and
            position of the solenoid relative to the origin of the world frame.
        """
        self.n = n
        self.I = I
        self.a = a
        self.b = b/2 # b parameter of model is half the length of the solenoid
        self.transform = transform

    def __call__(self, position, t):
        # Transform evaluation points so as to move solenoid to origin
        transformedPosition = self.transform.reverse(position)
        field = self.cartesianSolenoidModel(transformedPosition)
        # Rotate field vectors according to solenoid transform matrix
        return self.transform.apply(field)

    def cartesianSolenoidModel(self, position):
        """
        Compute field values of ideal solenoid given Cartesian co-ordinates.

        @param position: 3xN L{np.ndarray} of position co-ordinates.
        """
        x, y, z = position
        r,phi = cartesianToPolar(x,y)
        Br,Bz = self.cylindricalSolenoidModel(r,phi,z)
        Bx,By = polarToCartesian(Br,phi)
        return np.vstack((Bx,By,Bz))

    def cylindricalSolenoidModel(self,r,phi,z):
        """
        Compute field values of ideal solenoid given cylindrical co-ordinates.

        @param r: Radius value(s)
        @param phi: Angle value(s).
        @param z: Height value(s).
        """

        B_0 = 4*10**-7 * self.n * self.I
        z_plus = z+self.b
        z_minus = z-self.b
        z_plus2  = z_plus**2
        z_minus2 = z_minus**2

        alpha_plus = self.a/np.sqrt(z_plus2 + (r+self.a)**2)
        alpha_minus = self.a/np.sqrt(z_minus2 + (r+self.a)**2)
        beta_plus = z_plus/np.sqrt(z_plus2 + (r+self.a)**2)
        beta_minus = z_minus/np.sqrt(z_minus2 + (r+self.a)**2)
        gamma = (self.a-r)/(self.a+r)
        k_plus = np.sqrt((z_plus2 + (self.a-r)**2)/(z_plus2+(self.a+r)**2))
        k_minus = np.sqrt((z_minus2 + (self.a-r)**2)/(z_minus2+(self.a+r)**2))

        Br = B_0*(alpha_plus * _cel(k_plus,1,1,-1) -
                alpha_minus*_cel(k_minus,1,1,-1))
        Bz = ((B_0*self.a)/(self.a+r)) * \
                (beta_plus*_cel(k_plus,gamma**2,1,gamma) -
                beta_minus*_cel(k_minus,gamma**2,1,gamma))

        return Br,Bz

@np.vectorize
def _cel(kc,p,c,s):
    """
    Algorithm for the generalised complete elliptic integral C(k_c,p,c,s).
    """

    ERRTOL = 1e-6

    if kc == 0: return np.nan

    k = np.abs(kc)
    pp = p
    cc = c
    ss = s
    em = 1
    if p > 0:
        pp = np.sqrt(p)
        ss = s/pp
    else:
        f = kc**2
        q = 1-f
        g = 1-pp
        f = f-pp
        q = q*(ss-c*pp)
        pp = np.sqrt(f/g)
        cc = (c-ss)/g
        ss = -q/(g**2 * pp) + cc*pp

    f = cc
    cc = cc + ss/pp
    g = k/pp
    ss = 2*(ss+f*g)
    pp = g+pp
    g = em
    em = k+em
    kk = k

    while np.abs(g-k) > g*ERRTOL:
        k = 2*np.sqrt(kk)
        kk = k*em
        f = cc
        cc = cc + ss/pp
        g = kk/pp
        ss = 2*(ss + f*g)
        pp = g+pp
        g = em
        em = k+em

    return (np.pi/2)*(ss+cc*em)/(em*(em+pp))
