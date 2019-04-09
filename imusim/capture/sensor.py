"""
Classes to represent captured sensor data.
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
import pickle

class SensorDataCapture(object):
    """
    Sensor data captured synchronously from one or more devices.
    """

    def __init__(self):
        """
        Initialise capture.
        """
        self._devices = dict()

    def _addDevice(self, id, unit):
        self._devices[id] = unit

    @staticmethod
    def load(filename):
        """
        Load a capture from a file.

        @param filename: Name of file to load from.
        @return: A L{SensorDataCapture} object.
        """
        capture = pickle.load(open(filename, 'rb'), encoding='latin1')
        assert isinstance(capture, SensorDataCapture), \
                "File did not contain a SensorDataCapture"
        return capture

    def save(self,filename):
        """
        Save this capture to a file.

        @param filename: Name of file to save to.
        """
        pickle.dump(self, open(filename,'wb'))

    @property
    def devices(self):
        """
        List of the devices in this capture.
        """
        return self._devices.values()

    def device(self, id):
        """
        Get a device from the capture by identifier.

        @return: A L{CapturedDevice} object.
        """
        return self._devices[id]

class CapturedDevice(object):
    """
    A device with one or more colocated sensors from which data was captured.
    """
    def __init__(self, capture, id):
        """
        Initialise captured device.

        @param capture: Parent L{SensorDataCapture}.
        @param id: Device identifier.
        """
        self.capture = capture
        self.id = id
        self.capture._addDevice(id, self)
        self._sensorData = {}

    def addSensorData(self, id, data):
        """
        Add data for a sensor on this device.

        @param id: Identifier for this sensor.
        @param data: L{TimeSeries} of sensor data.
        """
        self._sensorData[id] = data

    @property
    def sensors(self):
        """
        List of the sensor identifiers on this sensor unit.
        """
        return self._sensorData.keys()

    def sensorData(self, id):
        """
        Get data for a given sensor identifier.

        @return: A L{TimeSeries} containing data for this sensor.
        """
        return self._sensorData[id]
