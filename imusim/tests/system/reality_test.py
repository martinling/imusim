"""
Test simulated outputs against real captured sensor data.
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
from imusim.io.qualisys_tsv import loadQualisysTSVFile
from imusim.capture.marker import SplinedMarkerCapture
from imusim.trajectories.multi_marker import MultiMarkerTrajectory
from imusim.trajectories.offset import OffsetTrajectory
from imusim.capture.sensor import SensorDataCapture
from imusim.trajectories.rigid_body import SampledBodyModel, SampledJoint
from imusim.trajectories.rigid_body import SplinedBodyModel
from imusim.platforms.imus import IdealIMU
from imusim.behaviours.imu import BasicIMUBehaviour
from imusim.environment.base import Environment
from imusim.maths.vector_fields import NaturalNeighbourInterpolatedField
from imusim.utilities.time_series import TimeSeries
from imusim.simulation.base import Simulation
from imusim.testing.vectors import assert_vectors_correlated
import numpy as np
from os import path

def testAgainstReality():
    dir = path.dirname(__file__)
    filebase = path.join(dir, "swing")
    refbase = path.join(dir, "stand")
    magbases = [path.join(dir, f) for f in ['magsweep1', 'magsweep2']]
    maglookup = {
        'Upper Leg IMU' : '66',
        'Orient 8' : '8',
        'Orient 43': '43'}
    magSamples = 2000
    refTime = 1.0
    posStdDev = 0.0005
    rotStdDev = 0.004
    ref3D = SplinedMarkerCapture(
        loadQualisysTSVFile(refbase + "_3D.tsv"), positionStdDev=posStdDev)
    ref6D = SplinedMarkerCapture(
        loadQualisysTSVFile(refbase + "_6D.tsv"), rotationStdDev=rotStdDev)
    capture3D = SplinedMarkerCapture(
        loadQualisysTSVFile(filebase + "_3D.tsv"), positionStdDev=posStdDev)
    captureSD = SensorDataCapture.load(filebase + ".sdc")
    hip, thigh, knee, shin, ankle = \
            ['Hip', 'Thigh', 'Knee Hinge', 'Shin', 'Ankle']
    jointNames = ['Upper Leg', 'Lower Leg', 'Foot']
    jointAbbrevs = ['femur', 'tibia', 'foot']
    orientIDs = ['66', '43', '8']
    jointMarkerNames = [hip, knee, ankle]
    refMarkerNames = [[thigh, knee], [shin, ankle], []]
    imuMarkerNames = \
            [[j + ' IMU - ' + str(i) for i in range(1,4)] for j in jointNames]
    jointMarkerSets = lambda c: [
        list(map(c.marker, jointMarkerNames)),
        [list(map(c.marker, r)) for r in refMarkerNames],
        [list(map(c.marker, i)) for i in imuMarkerNames]]
    imuMarkerSets = lambda c: [
        [c.marker(i[0]) for i in imuMarkerNames],
        [list(map(c.marker,i[1:])) for i in imuMarkerNames]]
    jointRefTrajectories = [MultiMarkerTrajectory(j, r + i, refTime=refTime)
        for j, r, i in zip(*(jointMarkerSets(ref3D)))]
    jointTrajectories = [
        MultiMarkerTrajectory(j, r + i, refVectors=m.refVectors) \
            for j, r, i, m in \
                zip(*(jointMarkerSets(capture3D) + [jointRefTrajectories]))]
    imuRefTrajectories = [MultiMarkerTrajectory(p, r, refTime=refTime)
        for p, r in zip(*(imuMarkerSets(ref3D)))]
    imuVecTrajectories = [MultiMarkerTrajectory(p, r, refVectors=m.refVectors)
        for p, r, m in zip(*(imuMarkerSets(capture3D) + [imuRefTrajectories]))]
    imuRefMarkers = [ref6D.marker(j + ' IMU') for j in jointNames]
    imuOffsets = [i.position(refTime) - j.position(refTime)
        for i, j in zip(imuRefTrajectories, jointRefTrajectories)]
    imuRotations = [t.rotation(refTime).conjugate * m.rotation(refTime)
        for t, m in zip(imuRefTrajectories, imuRefMarkers)]
    imuTrajectories = [OffsetTrajectory(v, o, r)
        for v, o, r in zip(imuVecTrajectories, imuOffsets, imuRotations)]
    imuData = [captureSD.device(i) for i in orientIDs]
    joints = []
    for i in range(len(jointNames)):
        name = jointNames[i]
        traj = jointTrajectories[i]
        if i == 0:
            model = SampledBodyModel(name)
            model.positionKeyFrames = traj.posMarker.positionKeyFrames
            joint = model
        else:
            parent = joints[-1]
            refTraj = jointRefTrajectories[i]
            parentRefTraj = jointRefTrajectories[i - 1]
            pos = refTraj.position(refTime)
            parentPos = parentRefTraj.position(refTime)
            joint = SampledJoint(joints[-1],name, offset=(pos - parentPos))
        joint.rotationKeyFrames = traj.rotationKeyFrames
        joints.append(joint)
    model = SplinedBodyModel(model)
    joints = model.joints
    imuJointTrajectories = [OffsetTrajectory(j, o, r)
        for j, o, r in zip(joints, imuOffsets, imuRotations)]
    positionSets = []
    valueSets = []
    for magbase in magbases:
        orient = SensorDataCapture.load(magbase + '.sdc')
        optical = loadQualisysTSVFile(magbase + '_6D.tsv')
        for marker in optical.markers:
            device = orient.device(maglookup[marker.id])
            magData = device.sensorData('magnetometer').values
            positionSets.append(marker.positionKeyFrames.values)
            valueSets.append(
                    marker.rotationKeyFrames.values.rotateVector(magData))
    positions = np.hstack(positionSets)
    values = np.hstack(valueSets)
    valid = ~np.any(np.isnan(positions),axis=0) & ~np.any(np.isnan(values),axis=0)
    dev = values - np.median(values[:,valid],axis=1).reshape((3,1))
    step = np.shape(values[:,valid])[1] // magSamples
    posSamples = positions[:,valid][:,::step]
    valSamples = values[:,valid][:,::step]
    environment = Environment()
    environment.magneticField = \
            NaturalNeighbourInterpolatedField(posSamples, valSamples)
    sim = Simulation(environment=environment)
    sim.time = model.startTime
    distortIMUs = []
    dt = 1/capture3D.sampled.frameRate
    for traj in imuJointTrajectories:
        platform = IdealIMU(sim, traj)
        distortIMUs.append(BasicIMUBehaviour(platform, dt))
    sim.run(model.endTime)
    for imu in range(3):
        for sensorName in ['accelerometer', 'magnetometer', 'gyroscope']:
            sim = getattr(distortIMUs[imu].imu,sensorName).rawMeasurements
            true = imuData[imu].sensorData(sensorName)(sim.timestamps + model.startTime)
            yield assert_vectors_correlated, sim.values, true, 0.8
