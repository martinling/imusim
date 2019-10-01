"""
Loader for Vicon optical capture data in CSV format.
"""

from __future__ import division
from imusim.capture.marker import MarkerCapture, Marker3DOF
import numpy as np
import csv

def loadViconCSVFile(filename):
    """
    Load 3DOF marker data from a Vicon CSV file.

    @param filename: Name of the CSV file to load.

    @return: A L{MarkerCapture} object.
    """

    captureRate = 60.0
    # Will change based on what Dan says
    data = []
    timestamps = []
    jointNames = []
    # A few empty lists
    reader = csv.reader(filename)
    # Read in the file
    names = [str(name) for name in reader[2]]
    # Set the third row of the file as the labels
    row = next(reader)
    # tmp = [str(element) for element in row[2:]]
    # Read in data starting from 3rd row of row (6th row in csv)

    for i, name in enumerate(names):
        if i%3 == 2:
            idx = name.rfind(':')
            name = name[idx+1:]
            jointNames.append(name)
    # Create list of jointnames from label

    for row in reader:
        timestamps.append(float(row[0]))
        # Set timestamp as first element of row
        newData = [float(element) for element in row[2:]]
        newData = np.split(np.array(newData), len(jointNames))
        data.append(newData)

    timestamps = np.array(timestamps)
    data = np.array(data)/1000.0
    # Data is in mm

    capture = MarkerCapture()
    capture.frameTimes = timestamps
    capture.frameCount = len(timestamps)
    capture.frameRate = captureRate
    capture.framePeriod = 1.0/captureRate
    
    for i, name in enumerate(jointNames):
        marker = Marker3DOF(capture, name)
        for time, position in zip(capture.frameTimes, data[:,i,:]):
            marker.positionKeyFrames.add(time, position.reshape(3,1))

    return capture
