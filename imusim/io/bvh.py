"""
Reading and writing BVH body model movement files.
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

from __future__ import with_statement
from imusim.trajectories.rigid_body import SampledBodyModel, SampledJoint, \
        PointTrajectory
from imusim.maths.quaternions import Quaternion
from imusim.maths.transforms import convertCGtoNED, convertNEDtoCG
import numpy as np

# Conversion co-efficient to convert CMU skeletal dimensions to meters
# Taken from the CMU FAQ http://mocap.cs.cmu.edu/faqs.php
CMU_CONVERSION = (1.0/0.45)*2.54/100.0

# Co-efficient to convert values in m to cm
M_TO_CM_CONVERSION = 100
CM_TO_M_CONVERSION = 0.01

def loadBVHFile(filename, conversionFactor=1):
    """
    Load sampled body model data from a BVH file.

    @param filename: Name of file to load.
    @param conversionFactor: Scale factor to apply to position values, e.g.
        to convert meters to centimeters use conversionFactor=100.

    @return: A L{SampledBodyModel} object.
    """
    with open(filename,'r') as bvhFile:
        loader = BVHLoader(bvhFile,conversionFactor)
        loader._readHeader()
        loader._readMotionData()
    return loader.model

class BVHLoader(object):
    '''
    Loader class for BVH files.
    '''
    def __init__(self, bvhFile, conversionFactor):
        '''
        Construct a new BVH loader for the specified BVH file.
        '''
        self.bvhFile = bvhFile
        self.tokens = []
        self.line = 0
        self.totalChannels = 0
        self.startOfMotionData = None
        self.conversionFactor = conversionFactor

    def _readHeader(self):
        '''
        Read the header information from the BVH file.
        '''
        assert self.bvhFile.tell() == 0, \
            'Must be at the start of file to read header'
        self._checkToken("HIERARCHY")
        self._checkToken("ROOT")
        self.model = SampledBodyModel()
        self.model.channels = []
        self.jointStack = []
        self.jointStack.append(self.model)
        self._readJoint()
        del self.jointStack

        self._checkToken("MOTION")
        self._checkToken("Frames:")
        self.frameCount = self._intToken()
        self._checkToken("Frame")
        self._checkToken("Time:")
        self.framePeriod = self._floatToken()

    def _readMotionData(self):
        '''
        Read the motion data from the BVH file.
        '''
        # update joint tree with new data.
        frame = 0
        lastRotation = {}
        for frame in range(self.frameCount):
            time = frame * self.framePeriod
            channelData = self._readFrame()
            for joint in self.model.joints:
                for chan in joint.channels:
                    chanData = channelData.pop(0)
                    if chan == "Xposition":
                        xPos = chanData
                    elif chan == "Yposition":
                        yPos = chanData
                    elif chan == "Zposition":
                        zPos = chanData
                    elif chan == "Xrotation":
                        xRot = chanData
                    elif chan == "Yrotation":
                        yRot = chanData
                    elif chan == "Zrotation":
                        zRot = chanData
                    else:
                        raise RuntimeError('Unknown channel: '+chan)

                if not joint.hasParent:
                    try:
                        p = np.array([[xPos,yPos,zPos]]).T * \
                                self.conversionFactor
                        p = convertCGtoNED(p)
                        del xPos, yPos, zPos
                        joint.positionKeyFrames.add(time, p)
                    except UnboundLocalError:
                        raise SyntaxError("No position data for root joint")

                try:
                    if joint.hasChildren:
                        q = Quaternion.fromEuler((zRot,xRot,yRot), order='zxy')
                        q = convertCGtoNED(q)
                        del xRot, yRot, zRot
                        # Rotation in BVH is relative to parent joint so
                        # combine the rotations
                        if joint.hasParent:
                            q = lastRotation[joint.parent] * q
                        joint.rotationKeyFrames.add(time, q)
                        lastRotation[joint] = q
                except UnboundLocalError:
                    raise SyntaxError("No rotation data for a joint that "
                            "is not an end effector")

    def _readFrame(self):
        '''
        Read a frame of data from the BVH file.

        @return: The channel data as a list of floats.
        '''
        chanData = self._readline().split()
        if len(chanData) != self.totalChannels:
            raise SyntaxError("Syntax error at line %d: Number of entries \
does not match the number of channels. (Found %d of %d)" %(
                self.line,len(chanData),self.totalChannels))
        return list(map(float,chanData))

    def _readJoint(self):
        # use jointStack[-1] to peek at the top of the jointStack
        currentJoint = self.jointStack[-1]
        name = self._token()
        # name will be 'Site' for End Site joints.
        if name == 'Site':
            name = currentJoint.parent.name+"_end"
        currentJoint.name = name

        self._checkToken("{")
        while True:
            token = self._token()
            if token == "OFFSET":
                x = self._floatToken()
                y = self._floatToken()
                z = self._floatToken()
                currentJoint.positionOffset = \
                        convertCGtoNED(np.array([[x,y,z]]).T*self.conversionFactor)
            elif token == "CHANNELS":
                n = self._intToken()
                channels = []
                for i in range(n):
                    token = self._token()
                    if token not in ["Xposition","Yposition","Zposition",\
                                     "Xrotation","Yrotation","Zrotation"]:
                        raise SyntaxError("Syntax error in line %d: Invalid \
channel name '%s'" %(self.line,token))
                    else:
                        channels.append(token)
                self.totalChannels += n
                currentJoint.channels = channels
            elif token in ("JOINT", "End"):
                if token == "JOINT":
                    joint = SampledJoint(currentJoint)
                else:
                    joint = PointTrajectory(currentJoint)
                joint.channels = []
                self.jointStack.append(joint)
                self._readJoint()
            elif token == "}":
                self.jointStack.pop()
                break
            else:
                raise SyntaxError("Syntax error in line %d: Unknown keyword '%s'" %(
                            self.line,token))

    def _checkToken(self,expectedToken):
        '''
        Try to read an expected token from the input queue.

        @raise SyntaxError: if the expected token is not found.
        '''
        token = self._token()
        if token != expectedToken:
            raise SyntaxError("Syntax error in line %d: Expected %s \
but found %s" %(self.line,expectedToken,token))

    def _intToken(self):
       '''
       Read an integer value from the token queue.
       '''
       token = self._token()
       try:
           return int(token)
       except ValueError:
           raise SyntaxError('Syntax error in line %d: Integer \
expected but found %s' %(self.line,token))

    def _floatToken(self):
        '''
        Read a float value from the token queue.
        '''
        token = self._token()
        try:
            return float(token)
        except ValueError:
            raise SyntaxError('Syntax error in line %d: Float \
expected but found %s' %(self.line,token))

    def _token(self):
        '''
        Return the next token in the input file.
        '''
        try:
            return self.tokens.pop(0)
        except IndexError:
            # There are no more tokens in the queue so read the next line
            l = self._readline()
            tokens = l.strip().split()
            for token in tokens:
                self.tokens.append(token)
            # recurse to return the next token or read another line
            return self._token()

    def _readline(self):
        '''
        Read a line from the input file and add the tokens to the list of
        tokens to process.
        '''
        l = self.bvhFile.readline()
        self.line+=1
        if l != "":
            return l
        else:
            raise StopIteration

def saveBVHFile(model, filename, samplePeriod, conversionFactor = 1):
    """
    Save a body model data to a BVH file.

    @param model: Body model to save.
    @param filename: Name of file to write to.
    @param samplePeriod: The time between frames in the BVH output.
    @param conversionFactor: Scale factor to apply to position values, e.g.
        to convert meters to centimeters use conversionFactor=100.
    """
    with open(filename,'w') as bvhFile:
        exporter = BVHExporter(model,bvhFile, samplePeriod, conversionFactor)
        exporter.writeHeader()
        for frame in range(exporter.frames):
            exporter.writeFrame(model.startTime + frame * samplePeriod)

class BVHExporter(object):
    """
    Exporter class to write BVH files.

    @ivar model: The body model to export.
    @ivar bvhFile: The BVH file to export data to.
    @ivar samplePeriod: The time between frames in the BVH output
    @ivar conversionFactor: Scale factor to apply to position values, e.g.
        to convert meters to centimeters use conversionFactor=100.
    """

    def __init__(self, model, bvhFile, samplePeriod, conversionFactor):
        '''
        Construct BVH exporter.

        @param model: Body model to save.
        @param bvhFile: File to write to.
        @param samplePeriod: The time between frames in the BVH output.
        @param conversionFactor: Scale factor to apply to position values, e.g.
            to convert meters to centimeters use conversionFactor=100.
        '''
        self.bvhFile = bvhFile
        self.model = model
        self.samplePeriod = samplePeriod
        self.conversionFactor = conversionFactor

    @property
    def frames(self):
        """ The number of frames that will be exported. """
        return int((self.model.endTime - self.model.startTime) //
                self.samplePeriod)

    def writeHeader(self):
        """
        Write the header of the BVH file.
        """
        self.bvhFile.write('HIERARCHY\n')
        pad = "  "
        def post(j):
            self.bvhFile.write("%s}\n" %(pad*j.depth))
        for p in self.model.preorderTraversal(postFunc=post):
            d = p.depth
            if not p.hasParent:
                self.bvhFile.write("%sROOT %s\n" %(pad*d, p.name))
            elif p.isJoint:
                self.bvhFile.write("%sJOINT %s\n" %(pad*d, p.name))
            else:
                self.bvhFile.write("%sEnd Site\n" %(pad*d))
            self.bvhFile.write("%s{\n" %(pad*d))

            self.bvhFile.write("%sOFFSET" %(pad*(d+1)))
            for v in convertNEDtoCG(p.positionOffset):
                self.bvhFile.write(" %.8f" % (v * self.conversionFactor))
            self.bvhFile.write("\n")

            if p.isJoint:
                if not hasattr(p, "channels"):
                    if not p.hasParent:
                        p.channels = ["Xposition", "Yposition", "Zposition",
                            "Zrotation", "Xrotation", "Yrotation"]
                    else:
                        p.channels = ["Zrotation", "Xrotation", "Yrotation"]

                self.bvhFile.write("%sCHANNELS %d" %(pad*(d+1),len(p.channels)))
                for chan in p.channels:
                    self.bvhFile.write(" %s" % chan)
                self.bvhFile.write("\n")

        self.bvhFile.write("MOTION\nFrames: %d\nFrame Time: %.8f\n" %(
            self.frames,self.samplePeriod))

    def writeFrame(self, t):
        for j in self.model.joints:
            if hasattr(j,'channels') and len(j.channels) > 0:

                if not j.hasParent:
                    pos = j.position(t) * self.conversionFactor
                    pos = convertNEDtoCG(pos)

                rot = j.rotation(t)
                if j.hasParent:
                    # BVH rotations are relative to parent joint
                    # RigidBodyModel rotations are absolute
                    rot = j.parent.rotation(t).conjugate * rot
                rot = convertNEDtoCG(rot)
                rot = rot.toEuler(order='zxy')

                for chan in j.channels:
                    if chan == 'Xposition':
                        self.bvhFile.write("%f " %pos[0])
                    elif chan == 'Yposition':
                        self.bvhFile.write("%f " %pos[1])
                    elif chan == 'Zposition':
                        self.bvhFile.write("%f " %pos[2])
                    elif chan == 'Zrotation':
                        self.bvhFile.write("%f " %rot[0])
                    elif chan == 'Xrotation':
                        self.bvhFile.write("%f " %rot[1])
                    elif chan == 'Yrotation':
                        self.bvhFile.write("%f " %rot[2])
        self.bvhFile.write("\n")
