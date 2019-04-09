"""
Reading of ASF/AMC body model movement files.
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

from imusim.trajectories.rigid_body import SampledBodyModel, SampledJoint
from imusim.trajectories.rigid_body import PointTrajectory
from imusim.maths.vectors import vector
from imusim.maths.transforms import convertCGtoNED
from imusim.maths.quaternions import Quaternion
from imusim.utilities.trees import TreeNode
from pyparsing import Word, Literal, Keyword, Group, Dict, Combine
from pyparsing import ZeroOrMore, OneOrMore, Suppress, Optional, LineEnd
from pyparsing import alphanums, nums, alphas, Regex, SkipTo

__all__ = ['loadASFFile']

class ASFBone(TreeNode):
    """
    Simple tree node for representing a bone in the ASF hierarchy.
    """
    def __init__(self, parent, bonedata, scale):
        """
        Consruct a new ASFBone.

        @param parent: The parent bone of this bone.
        @param bonedata: The L{ParseResults} object representing the bone
            data from the ASF file.
        """
        TreeNode.__init__(self, parent)
        self.name = bonedata.name
        axis = vector(*bonedata.axis)[::-1]
        rotOrder = bonedata.axisRotationOrder[::-1]
        self.rotationOffset = Quaternion.fromEuler(axis, rotOrder).conjugate
        self.childoffset = self.rotationOffset.rotateVector(
                vector(*bonedata.direction) * bonedata.length * scale)
        self.isDummy = bonedata.channels is ''
        self.channels = bonedata.channels

class ASFRoot(ASFBone):
    """
    Represents the root of an ASF model hierarchy
    """
    def __init__(self, rootdata):
        """
        Construct a new ASFRoot.

        @param rootdata: The L{ParseResults} object representing the root
            segment data from the ASF file.
        """
        TreeNode.__init__(self, None)
        self.name = 'root'
        self.childoffset = vector(0,0,0)
        self.isDummy = False
        axis = vector(*rootdata.axis)[::-1]
        rotOrder = rootdata.axisRotationOrder[::-1]
        self.rotationOffset = Quaternion.fromEuler(axis, rotOrder).conjugate
        self.channels = rootdata.channels

    def getBone(self, bonename):
        """
        Get a from the model tree by name.

        @param bonename: The name of the bone to return
        @return: The L{ASFBone} with the corresponding name, or `None` if such
            a bone does not exist.
        """
        for bone in self:
            if bone.name == bonename:
                return bone

def loadASFFile(asfFileName, amcFileName, scaleFactor, framePeriod):
    """
    Load motion capture data from an ASF and AMC file pair.

    @param asfFileName: Name of the ASF file containing the description of
        the rigid body model
    @param amcFileName: Name of the AMC file containing the motion data
    @param scaleFactor: Scaling factor to convert lengths to m. For data
        from the CMU motion capture corpus this should be 2.54/100 to
        convert from inches to metres.

    @return: A {SampledBodyModel} representing the root of the rigid body
        model structure.
    """
    with open(asfFileName, 'r') as asfFile:
        data = asfParser.parseFile(asfFile)
        scale = (1.0/data.units.get('length',1)) * scaleFactor

        bones = dict((bone.name,bone) for bone in data.bones)
        asfModel = ASFRoot(data.root)

        for entry in data.hierarchy:
            parent = asfModel.getBone(entry.parent)
            for childName in entry.children:
                ASFBone(parent, bones[childName], scale)

        imusimModel = SampledBodyModel('root')
        for subtree in asfModel.children:
            for bone in subtree:
                if not bone.isDummy:
                    offset = vector(0,0,0)
                    ancestors = bone.ascendTree()
                    while True:
                        ancestor = next(ancestors).parent
                        offset += ancestor.childoffset
                        if not ancestor.isDummy:
                            break

                    SampledJoint(parent=imusimModel.getJoint(ancestor.name),
                            name=bone.name,
                            offset=offset)
                if not bone.hasChildren:
                    PointTrajectory(
                            parent=imusimModel.getJoint(bone.name),
                            name=bone.name+'_end',
                            offset=bone.childoffset
                            )

        with open(amcFileName) as amcFile:
            motion = amcParser.parseFile(amcFile)
            t = 0
            for frame in motion.frames:
                for bone in frame.bones:
                    bonedata = asfModel.getBone(bone.name)
                    if bone.name == 'root':
                        data = dict((chan.lower(), v) for chan,v in
                            zip(bonedata.channels,bone.channels))
                        position = convertCGtoNED(scale * vector(data['tx'],
                            data['ty'],data['tz']))
                        imusimModel.positionKeyFrames.add(t, position)

                    axes, angles = zip(*[(chan[-1], angle) for chan, angle in
                        zip(bonedata.channels, bone.channels) if
                                chan.lower().startswith('r')])
                    rotation = (bonedata.rotationOffset.conjugate *
                            Quaternion.fromEuler(angles[::-1], axes[::-1]))
                    joint = imusimModel.getJoint(bone.name)
                    if joint.hasParent:
                        parentRot = joint.parent.rotationKeyFrames.latestValue
                        parentRotOffset = bonedata.parent.rotationOffset
                        rotation = parentRot * parentRotOffset * rotation
                    else:
                        rotation = convertCGtoNED(rotation)
                    joint.rotationKeyFrames.add(t, rotation)
                t += framePeriod

        return imusimModel

# Define parser tokens
comments = ZeroOrMore(Suppress(Literal('#') + SkipTo(LineEnd())))
intValue = Word(nums).setParseAction( lambda s,l,t: int(t[0]) )
floatValue = Regex(r'-?\d+(\.\d*)?(e-?\d*)?').setParseAction(lambda s,l,t: float(t[0]))
floatVector = Group(floatValue + floatValue + floatValue)
limit = Group(
        Suppress(Literal("(")) +
        floatValue +
        floatValue +
        Suppress(Literal(")")))
limits = Group(OneOrMore(limit))
channel = Word("TRtr","XYZxyz")
channels = Group(OneOrMore(channel))
rotationOrder = Word("XYZ", exact=3)
begin = Suppress(Keyword("begin"))
end = Suppress(Keyword("end"))
bonename = Combine(~end + Word(alphanums+"_-")).setWhitespaceChars(' ')

version = Keyword(":version") + Literal("1.10")
skeletonName = Keyword(":name") + bonename.setResultsName('name')
unitDefinition = Group(Word(alphas) + (floatValue | intValue | Word(alphas)))
unitSection = Keyword(":units") + \
        Dict(ZeroOrMore(unitDefinition)).setResultsName('units')
documentationSection = Keyword(':documentation') + \
        SkipTo(":").setResultsName('documentation')
rootSection = Group(Keyword(":root") &
        (Keyword("order") +
        channels.setResultsName('channels')) &
        (Keyword("position") +
        floatVector.setResultsName('position')) &
        (Keyword("axis") +
        rotationOrder.setResultsName("axisRotationOrder")) &
        (Keyword("orientation") +
        floatVector.setResultsName("axis"))
        ).setResultsName('root')
bone = Group(
        begin +
        Keyword("id") +
        intValue +
        Keyword("name") +
        bonename.setResultsName("name") +
        Keyword("direction") +
        floatVector.setResultsName("direction") +
        Keyword("length") +
        floatValue.setResultsName("length") +
        Keyword("axis") +
        floatVector.setResultsName("axis") +
        rotationOrder.setResultsName("axisRotationOrder") +
        Optional(
            Keyword("dof") +
            channels.setResultsName("channels") +
            Keyword("limits") +
            limits.setResultsName("limits")
            ) +
        end
        )

bonedataSection = (
        Keyword(":bonedata") +
        Group(ZeroOrMore(bone)).setResultsName("bones")
        )
hierarchyEntry = Group(bonename.setResultsName("parent") +
        Group(OneOrMore(bonename)).setResultsName("children")+
        Suppress(LineEnd()))
hierarchySection = (
        Keyword(":hierarchy") +
        begin + LineEnd() +
        Dict(OneOrMore(hierarchyEntry)).setResultsName("hierarchy") +
        end
        )

asfParser = (comments +
        version +
        skeletonName +
        unitSection +
        documentationSection +
        rootSection +
        bonedataSection +
        hierarchySection
        )

amcFrame = Group(intValue.setResultsName('frameNumber') + Suppress(LineEnd()) +
        OneOrMore(Group(
            bonename.setResultsName('name') +
            Group(OneOrMore(floatValue.setWhitespaceChars(' '))
                ).setResultsName('channels') +
            Suppress(LineEnd()))).setResultsName('bones'))
amcParser = (comments + SkipTo(intValue) +
        OneOrMore(amcFrame).setResultsName('frames'))
