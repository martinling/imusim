"""
Loader for Qualisys optical capture data in TSV format.
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
from imusim.capture.marker import MarkerCapture, Marker3DOF, Marker6DOF
from imusim.maths.quaternions import Quaternion, QuaternionArray
import numpy as np

def loadQualisysTSVFile(filename):
    """
    Load 3DOF or 6DOF marker data from a Qualisys Track Manager TSV file.

    @param filename: Name of TSV file to load.

    @return: A L{MarkerCapture} instance.
    """

    # Open file to read header.
    datafile = open(filename, 'r')

    # Count of header lines in file.
    headerLines = 0

    capture = MarkerCapture()

    markerIndices = dict()

    while True:

        # Read and count each file header line.
        line = datafile.readline().strip('\r\n')
        headerLines = headerLines + 1

        # Key and values are tab-separated.
        items = line.split('\t')
        key = items[0]
        values = items[1:]

        # Handle relevant fields.
        if key == 'FREQUENCY':
            # Frame rate.
            capture.frameRate = int(values[0])
            capture.framePeriod = 1/capture.frameRate
        elif key == 'DATA_INCLUDED':
            # 3D or 6D data.
            type = values[0]
            if (type == '3D'):
                # 3D data fields are implicitly XYZ co-ordinates.
                fieldNames = ['X', 'Y', 'Z']
        elif key == 'BODY_NAMES' or key == 'MARKER_NAMES':
            # List of markers or 6DOF body names.
            for i in range(len(values)):
                name = values[i]
                if type == '6D':
                    markerClass = Marker6DOF
                else:
                    markerClass = Marker3DOF
                markerIndices[markerClass(capture, name)] = i

            if (type == '6D'):
                # Field names are on a line after the body names, each
                # beginning with a space.
                fieldNames = [x.strip(' ') for x in
                    datafile.readline().strip('\r\n').split('\t')]
                headerLines = headerLines + 1

            # The above keys are always on the last line in the header.
            datafile.close()
            break

    # Load values from file, skipping header.
    data = np.loadtxt(filename, skiprows = headerLines)
    capture.frameCount = len(data)
    capture.frameTimes = np.arange(0, capture.frameCount / capture.frameRate,
            capture.framePeriod)

    # Find index in data array for a given marker and marker field name.
    def index(marker, name):
        offset = markerIndices[marker] * len(fieldNames)
        return offset + fieldNames.index(name)

    # Return data for a given marker and marker field names.
    def fields(marker, names):
        indices = [index(marker, name) for name in names]
        return data[:,indices[0]:indices[-1]+1]

    # Iterate through markers to fill in their data.
    for marker in capture.markers:
        # Get position data.
        positions = fields(marker, ['X', 'Y', 'Z']) / 1000
        if type == '6D':
            # Get rotation matrices.
            matrix_coeffs = fields(marker, ['Rot[%d]' % i for i in range(9)])
            # Mark invalid data (residual = -1) with NaNs.
            invalid = fields(marker, ['Residual'])[:,0] == -1.0
            matrix_coeffs[invalid] = np.nan
            positions[invalid] = np.nan
            # Convert to quaternions and unflip.
            quats = QuaternionArray(np.empty((capture.frameCount,4)))
            for i in range(len(quats)):
                q = Quaternion()
                q.setFromMatrix(matrix_coeffs[i].reshape((3,3)).T)
                quats[i] = q
            rotations = quats.unflipped()
            # Add rotation data to marker.
            for time, rotation in zip(capture.frameTimes, rotations):
                marker.rotationKeyFrames.add(time, rotation)
        # Add position data to marker.
        for time, position in zip(capture.frameTimes, positions):
            marker.positionKeyFrames.add(time, position.reshape(3,1))

    return capture
