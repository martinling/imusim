"""
Overall model of the environment experienced by simulated IMUs.
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

from imusim.maths.vector_fields import VectorField
from imusim.environment.magnetic_fields import EarthMagneticField
from imusim.environment.gravity import ConstantGravitationalField
from imusim.environment.radio_environment import RadioEnvironment, \
        IdealRadioEnvironment

class Environment(object):
    """
    Overall model of the environment experienced by simulated IMUs.

    An environment includes gravitational and magnetic vector fields, and a
    radio environment that models the propagation of radio transmissions.

    The default environment has constant magnetic and gravitational field
    models with nominal values for Edinburgh, UK, and an ideal, lossless
    radio environment.

    @ivar magneticField: L{VectorField} model of magnetic field.
    @ivar gravitationalField: L{VectorField} model of gravity.
    @ivar radioEnvironment: L{RadioEnvironment} model instance.
    """

    def __init__(self, magneticField=EarthMagneticField(),
            gravitationalField=ConstantGravitationalField(),
            radioEnvironment=IdealRadioEnvironment()):
        """
        Construct environment model.

        @param magneticField: L{VectorField} model of magnetic field.
        @param gravitationalField: L{VectorField} model of gravity.
        @param radioEnvironment: L{RadioEnvironment} model instance.
        """

        self.gravitationalField = gravitationalField
        self.magneticField = magneticField
        self.radioEnvironment = radioEnvironment
