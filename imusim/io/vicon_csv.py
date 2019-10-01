"""
Loader for Vicon optical capture data in CSV format.
"""

from __future__ import division
from imusim.capture.marker import MarkerCapture, Marker3DOF
import numpy as np

def loadViconCSVFile(filename):
    """
    Load 3DOF marker data from a Vicon CSV file.

    @param filename: Name of the CSV file to load.

    @return: A L{MarkerCapture} object.
    """
    # Open file to read header.
    datafile = open(filename, 'r')

    # Function to get comma-separated values from a line.
    values = lambda line: line.rstrip('\r\n').split(',')
    # Read the word "Trajectories".
    datafile.readline()
    # Read the number of trajectories.
    trajectories = datafile.readline()
    # Get the unstripped marker names.
    markers = values(datafile.readline())
    # Strip the markers to only get every 3rd item starting from index 2.
    markernames = [n.split(':')[-2] for n in markers[2::3]]
    # Get the column names.
    colnames = datafile.readline()
    # Read the line containing the units of measurement (mm by default).
    datafile.readline()

    # Get data.
    data = np.array([[float(v or np.nan) for v in values(line)]
        for line in datafile.readlines()]).T

    frameTimes = data[1]
    positions = data[2:].reshape((len(markernames),3,-1)) / 1000

    capture = MarkerCapture()
    capture.frameTimes = frameTimes
    capture.frameCount = len(frameTimes)
    capture.framePeriod = np.min(np.diff(frameTimes))
    capture.frameRate = 1 / capture.framePeriod

    for i, name in enumerate(markernames):
        marker = Marker3DOF(capture, name)
        for time, position in zip(capture.frameTimes, positions[i].T):
            marker.positionKeyFrames.add(time, position.reshape(3,1))

    return capture
