"""
Base class for MAC algorithm implementations.
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
from imusim.platforms.base import Platform
from imusim.platforms.radios import Radio, RadioPacket
from imusim.behaviours.timing import TimerMultiplexer
from imusim.maths.quaternions import Quaternion
from imusim.utilities.documentation import prepend_method_doc
import numpy as np


class MAC(ABC):
    """
    Base class for MAC implementations.

    @ivar radio: The L{Radio} used by the MAC.
    @ivar timerMux: L{TimerMultiplexer} used by the MAC.
    """
    def __init__(self, radio, timerMux):
        """
        Initialise MAC.

        @param radio: The L{Radio} to be used by the MAC.
        @param timerMux: L{TimerMultiplexer} to be used by the MAC.
        """
        self.radio = radio
        self.timerMux = timerMux

        self.radio.setReceiveHandler(self.handleIncomingPacket)

    @abstractmethod
    def queuePacket(self, packet):
        """
        Queue a packet for transmission

        @param packet: The L{RadioPacket} to transmit.
        """
        pass

    @abstractmethod
    def handleIncomingPacket(self, packet):
        """
        Handle a received packet from the radio.

        @param packet: the received L{RadioPacket}.
        """
        pass
