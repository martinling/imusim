"""
Test of simulation with distributed linear acceleration estimator.
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

from imusim.simulation.base import Simulation
from imusim.io.bvh import loadBVHFile
from imusim.algorithms.orientation import DistLinAccelCF
from imusim.platforms.imus import IdealIMU
from imusim.simulation.calibrators import ScaleAndOffsetCalibrator
from imusim.behaviours.imu import BasicIMUBehaviour
from imusim.behaviours.timing import TimerMultiplexer
from imusim.behaviours.orient_tdma import InterSlaveSchedule
from imusim.behaviours.orient_tdma import MasterMAC
from imusim.behaviours.orient_tdma import InterSlaveMAC
from imusim.behaviours.orient_tdma import DataPacket
from imusim.behaviours.orient_tdma import AuxPacket
from imusim.trajectories.base import StaticTrajectory
from imusim.trajectories.offset import OffsetTrajectory
from imusim.platforms.basestations import IdealBasestation
from imusim.trajectories.rigid_body import SampledBodyModel
from imusim.trajectories.rigid_body import SplinedBodyModel
from imusim.behaviours.basestation import BodyModelReconstructor
from imusim.testing.quaternions import assert_quaternions_correlated
from imusim.testing.random_data import randomPosition
from imusim.utilities.time_series import TimeSeries
from os import path
import numpy as np

def testDistLinAccSimulation():
    sim = Simulation()
    samplingPeriod = 1.0/100
    calibrator = ScaleAndOffsetCalibrator(sim.environment, 1000,
            samplingPeriod, 20)
    bvhFile = path.join(path.dirname(__file__), 'walk.bvh')
    sampledModel = loadBVHFile(bvhFile, conversionFactor=0.01)
    posTimes = sampledModel.positionKeyFrames.timestamps
    sampledModel.positionKeyFrames = TimeSeries(
            posTimes, np.zeros((3,len(posTimes))))
    simModel = SplinedBodyModel(sampledModel)
    joints = list(simModel.joints)
    sim.time = simModel.startTime
    k = 128
    slotTime = 0.001
    txSlots = list(range(2, len(joints))) + [0,1]
    auxRxSlots = range(1, len(joints))
    auxTxSlots = [joints.index(j.parent) for j in joints[1:]]
    schedule = InterSlaveSchedule(slotTime, slotTime, txSlots,
            auxTxSlots, auxRxSlots)

    def setupIMU(id, joint):
        offset = randomPosition((-0.1, 0.1))
        platform = IdealIMU()
        calibration = calibrator.calibrate(platform)
        platform.simulation = sim
        platform.trajectory = OffsetTrajectory(joint, offset)
        filter = DistLinAccelCF(samplingPeriod,
                platform.trajectory.rotation(simModel.startTime),
                k, joint, offset)

        def updateChildren():
            for child in filter.children:
                childID = joints.index(child)
                childAccel = filter.childAcceleration(child.positionOffset, samplingPeriod)
                if childAccel is not None:
                    auxPacket = AuxPacket(id, childID,
                            [('linearAcceleration', childAccel)])
                    mac.queueAuxPacket(auxPacket)

        def handlePacket(packet):
            filter.handleLinearAcceleration(packet['linearAcceleration'], samplingPeriod)
            updateChildren()

        def handleSample(behaviour):
            rotation = filter.rotation.latestValue
            packet = DataPacket(id, [('rotation', rotation)])
            mac.queuePacket(packet)
            if not filter.joint.hasParent:
                updateChildren()

        behaviour = BasicIMUBehaviour(platform, samplingPeriod, calibration,
                filter, handleSample, initialTime=sim.time)
        mac = InterSlaveMAC(platform.radio, behaviour.timerMux, schedule,
                id, handlePacket)

    imus = [setupIMU(i, joint) for i, joint in enumerate(joints)]

    basestation = IdealBasestation(sim, StaticTrajectory())
    recModel = SampledBodyModel.structureCopy(simModel)
    reconstructor = BodyModelReconstructor(recModel, initialTime=sim.time)

    def handleFrame(packets):
        for packet in packets:
            packet['jointName'] = recModel.jointNames[packet.source]
        reconstructor.handleData(packets, schedule.framePeriod)

    basestation.radio.channel = schedule.dataChannel
    MasterMAC(basestation.radio, TimerMultiplexer(basestation.timer),
            schedule, handleFrame)

    sim.run(simModel.endTime)

    for simJoint, recJoint in zip(joints, recModel.joints):
        times = simJoint.rotationKeyFrames.timestamps
        assert_quaternions_correlated(simJoint.rotation(times),
            recJoint.rotation(times), 0.85)
