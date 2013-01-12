"""
Tests for overall system simulation.
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
from imusim.algorithms import orientation
from imusim.algorithms.orientation import OrientationFilter
from imusim.testing.inspection import getImplementations
from imusim.platforms.imus import MagicIMU
from imusim.simulation.calibrators import ScaleAndOffsetCalibrator
from imusim.behaviours.imu import BasicIMUBehaviour
from imusim.behaviours.timing import TimerMultiplexer
from imusim.behaviours.orient_tdma import MasterMAC, SlaveMAC
from imusim.behaviours.orient_tdma import Schedule, DataPacket
from imusim.trajectories.base import StaticTrajectory
from imusim.platforms.basestations import IdealBasestation
from imusim.trajectories.rigid_body import SampledBodyModel, SplinedBodyModel
from imusim.behaviours.basestation import BodyModelReconstructor
from imusim.testing.quaternions import assert_quaternions_correlated
from os import path
import numpy as np

testMotion = path.join(path.dirname(__file__), 'walk.bvh')

filterParameters = {
    orientation.OrientCF : {
        'k' : 1,
        'aT' : 1},
    orientation.BachmannCF : {
        'k' : 1},
    orientation.YunEKF : {
        'initialCovariance' : np.asmatrix(
            np.diag([0.01]*3 + [0.0001]*4)),
        'measurementCovariance' : np.asmatrix(
            np.diag([0.01]*3 + [0.0001]*4)),
        'D' : 50,
        'tau' : 0.5}
    }

def checkSystemSimulation(filterClass):
    sim = Simulation()
    env = sim.environment
    samplingPeriod = 1.0/100
    calibrator = ScaleAndOffsetCalibrator(env, 1000, samplingPeriod, 20)
    simModel = SplinedBodyModel(loadBVHFile(testMotion, conversionFactor=0.01))
    sim.time = simModel.startTime
    slotTime = 0.001
    schedule = Schedule(slotTime, slotTime, range(len(list(simModel.joints))))

    def setupIMU(id, joint):
        platform = MagicIMU()
        calibration = calibrator.calibrate(platform)
        platform.simulation = sim
        platform.trajectory = joint
        filter = filterClass(initialRotation=joint.rotation(simModel.startTime),
                initialTime=sim.time,
                initialRotationalVelocity=joint.rotation(simModel.startTime)
                    .rotateFrame(joint.rotationalVelocity(simModel.startTime)),
                gravFieldReference=env.gravitationalField.nominalValue,
                magFieldReference=env.magneticField.nominalValue,
                **filterParameters.get(filterClass, {}))
        def handleSample(behaviour):
            rotation = filter.rotation.latestValue
            packet = DataPacket(id, [('rotation', rotation)])
            behaviour.mac.queuePacket(packet)
        behaviour = BasicIMUBehaviour(platform, samplingPeriod, calibration,
                filter, handleSample, initialTime=sim.time)
        behaviour.mac = SlaveMAC(platform.radio, behaviour.timerMux,
                schedule, id)
        return behaviour

    imus = [setupIMU(i, joint) for i, joint in enumerate(simModel.joints)]

    basestation = IdealBasestation(sim, StaticTrajectory())
    recModel = SampledBodyModel.structureCopy(simModel)
    reconstructor = BodyModelReconstructor(recModel, initialTime=sim.time)

    def handleFrame(packets):
        for packet in packets:
            packet['jointName'] = recModel.jointNames[packet.source]
        reconstructor.handleData(packets, schedule.framePeriod)

    MasterMAC(basestation.radio, TimerMultiplexer(basestation.timer),
            schedule, handleFrame)
    sim.run(simModel.endTime)
    for simJoint, recJoint in zip(simModel.joints, recModel.joints):
        times = simJoint.rotationKeyFrames.timestamps
        assert_quaternions_correlated(simJoint.rotation(times),
            recJoint.rotation(times), 0.9)

def testSystemSimulation():
    for filterClass in getImplementations(orientation, OrientationFilter):
        if filterClass is not orientation.DistLinAccelCF:
            yield checkSystemSimulation, filterClass
