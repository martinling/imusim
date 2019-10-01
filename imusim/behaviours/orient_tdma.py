"""
Simulation model of the TDMA scheme used by the Orient IMU system.
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
from imusim.behaviours.mac import MAC
from imusim.behaviours.timing import VirtualTimer
from imusim.platforms.radios import RadioPacket
from imusim.maths.quaternions import Quaternion
from imusim.utilities.documentation import prepend_method_doc
import numpy as np

class Schedule:
    """
    Transmission schedule.
    """

    def __init__(self, dataSlotTime, syncSlotTime, dataSlots):
        """
        Initialise schedule.

        @param dataSlotTime: Data slot duration (float, s).
        @param syncSlotTime: Synchronisation slot duration (float, s).
        @param dataSlots: List of device IDs indexed by allocated slot.
        """
        self.dataSlotTime = dataSlotTime
        self.syncSlotTime = syncSlotTime
        self.dataSlots = dataSlots

    @property
    def framePeriod(self):
        """ Duration of one complete frame. """
        return len(self.dataSlots) * self.dataSlotTime + self.syncSlotTime

    def dataTxSlot(self, id):
        """ Transmission slot number for given device ID. """
        return self.dataSlots.index(id)

class SyncPacket(RadioPacket):
    """
    Synchronisation packet from master.
    """
    @property
    def bytes(self):
        return 1

class DataPacket(RadioPacket):
    """
    Data packet from slave to master.

    @ivar source: Source device ID.
    """
    def __init__(self, source, data):
        """
        Initialise data packet.
        @param source: Source device ID.
        @param data: Payload data.
        """
        RadioPacket.__init__(self, data)
        self.source = source

    @property
    def bytes(self):
        size = 2
        for key in self.keys():
            data = self[key]
            if isinstance(data, Quaternion):
                size += array(data.components).nbytes
            elif isinstance(data, numpy.ndarray):
                size += data.nbytes
        return size

class AuxPacket(DataPacket):
    """
    Auxillary packet from slave to slave.

    @ivar source: Source device ID.
    @ivar dest: Destination device ID.
    """
    def __init__(self, source, dest, data):
        """
        Initialise auxillary packet.

        @param source: Source device ID.
        @param dest: Destination device ID.
        @param data: Payload data.
        """
        DataPacket.__init__(self, source, data)
        self.dest = dest

class MasterMAC(MAC):
    """
    Master MAC implementation.
    """

    @prepend_method_doc(MAC)
    def __init__(self, radio, timerMux, schedule, frameHandler=None):
        """
        @param schedule: L{Schedule} object giving transmission scheduling.
        @param frameHandler: Handler function to call at the end of each frame.
            The packets received in that frame will be passed as a list.
        """
        MAC.__init__(self, radio, timerMux)
        self.schedule = schedule
        self.frameHandler = frameHandler
        self.frameTimer = VirtualTimer(timerMux, self._frameTimerHandler)
        self.packets = []
        self.frameTimer.start(self.schedule.framePeriod, repeat=True)

    def _frameTimerHandler(self):
        if self.frameHandler is not None:
            self.frameHandler(self.packets)
        self.packets = []
        self.radio.transmit(SyncPacket())

    def queuePacket(self, packet):
        raise NotImplementedError

    def handleIncomingPacket(self, packet):
        if isinstance(packet, DataPacket):
            self.packets.append(packet)

class SlaveMAC(MAC):
    """
    Slave MAC implementation.
    """

    @prepend_method_doc(MAC)
    def __init__(self, radio, timerMux, schedule, id):
        """
        @param schedule: L{Schedule} object giving transmission scheduling.
        @param id: The device ID of this node.
        """
        MAC.__init__(self, radio, timerMux)
        self.schedule = schedule
        self.id = id
        self.dataTxSlot = self.schedule.dataTxSlot(self.id)
        self.slotTimer = VirtualTimer(timerMux, self._slotTimerHandler)
        self.slot = len(self.schedule.dataSlots)
        self.txPacket = None

    def _slotTimerHandler(self):
        if self.slot == len(self.schedule.dataSlots):
            # End of sync slot, start of first data slot.
            self.slot = 0
        else:
            self.slot += 1

        if self.slot == len(self.schedule.dataSlots):
            # Start of sync slot.
            self.slotTimer.start(self.schedule.syncSlotTime)
        else:
            self.slotTimer.start(self.schedule.dataSlotTime)

        if self.slot == self.dataTxSlot:
            # Our transmit slot.
            if self.txPacket is not None:
                self.radio.transmit(self.txPacket)
                self.txPacket = None

    def queuePacket(self, packet):
        if isinstance(packet, DataPacket):
            self.txPacket = packet
        else:
            raise ValueError

    def handleIncomingPacket(self, packet):
        if isinstance(packet, SyncPacket):
            self.slotTimer.start(self.schedule.syncSlotTime/2)

class InterSlaveSchedule(Schedule):
    """
    Transmission schedule including an auxillary slave-to-slave channel.
    """

    @prepend_method_doc(Schedule)
    def __init__(self, dataSlots, dataSlotTime, syncSlotTime, auxTxSlots, auxRxSlots, dataChannel=0, auxChannel=1):
        """
        @param auxTxSlots: List of device IDs indexed by allocated
            auxillary TX slot.
        @param auxRxSlots: List of device IDs indexed by allocated
            auxillary RX slot.
        @param dataChannel: Channel identifier of master-slave channel.
        @param auxChannel: Channel identifier of slave-slave channel.
        """
        Schedule.__init__(self, dataSlots, dataSlotTime, syncSlotTime)
        self.auxTxSlots = auxTxSlots
        self.auxRxSlots = auxRxSlots
        self.dataChannel = dataChannel
        self.auxChannel = auxChannel

    def _auxSlots(self, sequence, id):
        return list(np.nonzero(np.array(sequence)==id)[0])

    def slaveAuxTxSlots(self, id):
        return self._auxSlots(self.auxTxSlots, id)

    def slaveAuxRxSlots(self, id):
        return self._auxSlots(self.auxRxSlots, id)

    def auxSlot(self, source, dest):
        return np.nonzero((np.array(self.auxTxSlots)==source)
                & (np.array(self.auxRxSlots)==dest))[0][0]

class InterSlaveMAC(SlaveMAC):
    """
    Slave MAC implementation with support for slave-slave packets.
    """

    @prepend_method_doc(SlaveMAC)
    def __init__(self, radio, timerMux, schedule, id, packetHandler=None):
        """
        @param packetHandler: Handler to call on reception of a slave-slave
            packet. The packet will be passed as an argument.
        """
        SlaveMAC.__init__(self, radio, timerMux, schedule, id)
        self.auxTxSlots = self.schedule.slaveAuxTxSlots(self.id)
        self.auxRxSlots = self.schedule.slaveAuxRxSlots(self.id)
        self.auxSlotTimer = VirtualTimer(timerMux, self._auxSlotTimerHandler)
        self.auxPackets = {}
        self.packetHandler = packetHandler
        self.radio.channel = schedule.dataChannel

    def _slotTimerHandler(self):
        SlaveMAC._slotTimerHandler(self)
        if self.slot == 0:
            self.auxSlot = 0
        else:
            self.auxSlot += 1
        if self.slot in self.schedule.dataSlots:
            self.auxSlotTimer.start(self.schedule.dataSlotTime/2)
            if self.auxSlot in self.auxTxSlots + self.auxRxSlots:
                self.radio.channel = self.schedule.auxChannel

    def _auxSlotTimerHandler(self):
        if self.auxSlot in self.auxTxSlots:
            packet = self.auxPackets.pop(self.auxSlot, None)
            if packet is not None:
                self.radio.transmit(packet)
        if self.slot + 1 in (self.dataTxSlot, len(self.schedule.dataSlots)):
            self.radio.channel = self.schedule.dataChannel

    def queueAuxPacket(self, packet):
        """
        Queue a packet for transmission to another slave.

        @param packet: L{RadioPacket} to transmit.
        """
        if isinstance(packet, AuxPacket):
            slot = self.schedule.auxSlot(packet.source, packet.dest)
            self.auxPackets[slot] = packet
        else:
            raise ValueError

    def handleIncomingPacket(self, packet):
        SlaveMAC.handleIncomingPacket(self, packet)
        if isinstance(packet, AuxPacket) and packet.dest == self.id:
            if self.packetHandler is not None:
                self.packetHandler(packet)
