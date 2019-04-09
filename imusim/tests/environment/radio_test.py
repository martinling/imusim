"""
Tests for radio environments.
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

from imusim.environment.radio_environment import BERRadioEnvironment
from imusim.platforms.radios import IdealRadio, RadioPacket
from imusim.platforms.base import Platform
from imusim.simulation.base import Simulation
from imusim.environment.base import Environment

class TestPlatform(Platform):
    def __init__(self, simulation=None, trajectory=None):
        self.radio = IdealRadio(self)
        self.packetsReceived = []
        self.radio.setReceiveHandler(self.handlePacket)
        Platform.__init__(self, simulation, trajectory)

    def handlePacket(self, packet):
        self.packetsReceived.append(packet)

    def sendPacket(self, packet):
        self.radio.transmit(packet)

    @property
    def components(self):
        return [self.radio]

class TestPacket(RadioPacket):
    @property
    def bytes(self):
        return 1

def testIdealRadioEnvironment():
    sim = Simulation()
    tx = TestPlatform(sim)
    rx = TestPlatform(sim)
    packet = TestPacket()

    tx.sendPacket(packet)

    assert len(tx.packetsReceived) == 0
    assert len(rx.packetsReceived) == 1

def testBERRadioEnvironment():
    env = Environment(radioEnvironment=BERRadioEnvironment(1e-4, seed=0))

    sim = Simulation(environment=env)
    tx = TestPlatform(sim)
    rx = TestPlatform(sim)
    packet = TestPacket()

    for _ in range(1000):
        tx.sendPacket(packet)

    assert len(tx.packetsReceived) == 0
    assert 0 < len(rx.packetsReceived) < 1000


