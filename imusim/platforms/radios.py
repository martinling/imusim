"""
Radio models.
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
from imusim.platforms.base import Component


class RadioPacket(dict, ABC):
    """
    Class to represent a radio packet.

    Payload fields should be added as dictionary items.

    Metadata relevant to radio simulation should be added as attributes.
    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        """
        Construct a radio packet.
        """
        pass

    @property
    @abstractmethod
    def bytes(self):
        pass


class Radio(Component):
    """
    Base class for radios.
    """
    def __init__(self, platform):
        Component.__init__(self, platform)
        self._receiveHandler = None

    def setReceiveHandler(self, receiveHandler):
        """
        Set the packet handler function.

        @param receiveHandler: Function to be called when a packet is
            successfully received. The received L{RadioPacket} will be passed
            as an argument.
        """
        self._receiveHandler = receiveHandler

    @abstractmethod
    def handlePacket(self, packet):
        """
        Simulate the result of an incoming packet from the radio channel. If
        the packet is received, it should be passed to the receive handler.
        """
        pass

    @abstractmethod
    def transmit(self, packet):
        """
        Transmit a packet.

        @param packet: L{RadioPacket} to send.
        """
        pass

    @property
    def _env(self):
        """
        Shortcut to access the L{RadioEnvironment} experienced by this radio.
        """
        return self.platform.simulation.environment.radioEnvironment

    def _simulationChange(self):
        if self.platform.simulation is not None:
            self._env.receivers.add(self)


class IdealRadio(Radio):
    """
    An ideal radio that can always receive and transmits instantaneously.
    """
    def handlePacket(self, packet):
        if self._receiveHandler is not None:
            self._receiveHandler(packet)

    def transmit(self, packet):
        self._env.transmit(self, packet)
