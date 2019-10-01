"""
Analogue to digital converter models.
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

from abc import abstractmethod
from imusim.simulation.base import Simulation
from imusim.platforms.base import Component
from imusim.platforms.sensors import Sensor
from imusim.utilities.documentation import prepend_method_doc
import SimPy.Simulation
import numpy as np

class ADC(Component):
    """
    Base class for ADCs.
    """
    @abstractmethod
    def transferFunction(self, voltageValues):
        """
        Apply the ADC transfer function that converts voltages to ADC points.

        @param voltageValues: 3xN L{np.ndarray} of input voltages.
        @return: 3xN L{np.ndarray} of quantised output values.
        """
        pass

    @abstractmethod
    def startSample(self, callback, *sensors):
        """
        Initiate sampling of sensors.

        @param callback: Function to call when sampling is complete. The
            function should accept a list of sampled values returned in the
            same order as the sensors where specified in this call.
        @param sensors: List of L{Sensor} objects to sample, in order.
        """
        pass

class IdealADC(ADC):
    """
    An ideal ADC that has no quantisation or delay between samples.
    """

    def transferFunction(self, voltageValues):
        return voltageValues

    def startSample(self, callback, *sensors):
        data = []
        for sensor in sensors:
            data.append(self.transferFunction(
                sensor.voltages(self.platform.simulation.time)))

        callback(*data)

class QuantisingADC(ADC):
    """
    A quantising ADC with configurable resolution.

    The ADC has a transfer function given by::

                   / -2^(bits-1),                       V < -Vref
      adc points = | floor(V/Vref * 2^(bits-1) + 0.5) - 2^(bits-1)
                   \\ 2^(bits-1) -1                     V >= Vref
    """

    def __init__(self, platform, bits, vref):
        """
        @param bits: Resolution of the ADC in output bits
        @param vref: Reference voltage used for voltage conversion
        """
        ADC.__init__(self, platform)
        self._bits = bits
        self._factor = 2**(bits-1)
        self._vref = vref

    def transferFunction(self, voltageValues):
        return np.clip(np.floor((voltageValues/self._vref) * self._factor +
            0.5) - self._factor, -self._factor , self._factor - 1)

    def startSample(self, callback, *sensors):
        data = []
        for sensor in sensors:
            data.append(self.transferFunction(
                sensor.voltages(self.platform.simulation.time)))

        callback(*data)

class ParametricADC(QuantisingADC):
    """
    A parametric ADC with configurable resolution and delay between samples.
    """

    @prepend_method_doc(QuantisingADC)
    def __init__(self, platform, bits, vref, sampleDelay):
        """
        @param sampleDelay: Time in seconds between consecutive samples
        """
        QuantisingADC.__init__(self, platform, bits, vref)
        self._sampleDelay = sampleDelay

    class _SampleProcess(SimPy.Simulation.Process):
        def __init__(self, adc, callback, *sensors):
            SimPy.Simulation.Process.__init__(self,
                    sim=adc.platform.simulation.engine)
            self.adc = adc
            self.callback = callback
            self.sensors = sensors
            self.sim.activate(self, self.eventLoop())

        def eventLoop(self):
            data = []
            for sensor in self.sensors:
                voltages = np.empty((3, 1))
                for axis in [0,1,2]:
                    yield SimPy.Simulation.hold, self, self.adc._sampleDelay
                    voltages[axis] = sensor.voltages(
                            self.adc.platform.simulation.time)[axis]
                data.append(self.adc.transferFunction(voltages))
            self.callback(*data)

    def startSample(self, callback, *sensors):
        self._SampleProcess(self, callback, *sensors)
