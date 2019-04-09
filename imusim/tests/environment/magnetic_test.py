"""
Tests for magnetic field environment.
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

from imusim.environment.magnetic_fields import EarthMagneticField
from imusim.environment.magnetic_fields import DistortedMagneticField
from imusim.environment.magnetic_fields import SolenoidMagneticField
from numpy.testing import assert_almost_equal
from imusim.maths.vectors import vector
from imusim.maths.transforms import AffineTransform
from numpy import random

# Tuple of constructor args and nominal field vectors for Earth fields
TEST_PARAMS = [
        ((1, 0, 0), vector(1,0,0)),
        ((1, 90, 0), vector(0,0,1)),
        ((1, 0, 90), vector(0,1,0))
        ]

def checkBaseField(field, nominal):
    assert_almost_equal(field.nominalValue, nominal)
    assert_almost_equal(field(vector(0,0,0), 0), nominal)

def testEarthField():
    for args, nominal in TEST_PARAMS:
        yield checkBaseField, EarthMagneticField(*args), nominal

def checkZeroHeadingVariation(field, position):
    assert_almost_equal(field.headingVariation(position, 0), 0)

def testUndistortedField():
    random.seed(0)
    for args, nominal in TEST_PARAMS:
        base = EarthMagneticField(*args)
        field = DistortedMagneticField(base)
        yield checkBaseField, field, nominal
        for i in range(10):
            yield checkZeroHeadingVariation, field, random.uniform(size=(3,1),
                    low=-10, high=10)

def checkNonZeroHeadingVariation(field, position):
    assert abs(field.headingVariation(position, 0)) > 0.01

def testDistortedField():
    random.seed(0)
    field = DistortedMagneticField(EarthMagneticField())
    field.addDistortion(SolenoidMagneticField(200,20,0.05,0.2,
        AffineTransform()))
    for i in range(10):
        yield checkNonZeroHeadingVariation, field, random.uniform(size=(3,1),
                low=-1, high=1)



