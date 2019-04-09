"""
Radio environment models.
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
from abc import ABC, abstractmethod


class RadioEnvironment(ABC):
    """
    Base class for radio environment models.

    @ivar receivers: list of potential receivers
    """
    def __init__(self):
        """
        Construct radio environment model.
        """
        self.receivers = set()

    @abstractmethod
    def transmit(self, transmitter, packet):
        """
        Simulate transmitting a packet in the environment.

        @param transmitter: The L{Radio} sending the packet.
        @param packet: The L{RadioPacket} being sent.
        """
        pass


class IdealRadioEnvironment(RadioEnvironment):
    """
    An ideal radio channel where all packet transmissions are successful.
    """
    def __init__(self):
        RadioEnvironment.__init__(self)

    def transmit(self, transmitter, packet):
        txChannel = getattr(transmitter, 'channel', None)
        for receiver in self.receivers:
            rxEnabled = getattr(receiver, '_receiverEnabled', True)
            rxChannel = getattr(receiver, 'channel', None)
            if rxEnabled and txChannel == rxChannel and receiver is not transmitter:
                receiver.handlePacket(packet)


class BERRadioEnvironment(IdealRadioEnvironment):
    """
    Radio channel with a constant bit error rate.

    Packets will be dropped at random with a probability dependent on their
    length and the BER of the environment.
    """
    def __init__(self, ber, seed=None):
        """
        Construct BER radio environment model.

        @param ber: Probability of error in a single bit.
        @param seed: Seed value for randomisation.
        """
        RadioEnvironment.__init__(self)
        self.ber = ber
        self.rng = np.random.RandomState(seed=seed)

    def transmit(self, transmitter, packet):
        lossProbability = 1 - (1 - self.ber)**(packet.bytes*8)
        if self.rng.uniform(0,1) > lossProbability:
            IdealRadioEnvironment.transmit(self, transmitter, packet)
