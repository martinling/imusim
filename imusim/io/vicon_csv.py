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

    data = []
    timestamps = []
    jointNames = []
    raw = []
    # A few empty lists
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        # Read in the file
        for row in reader:
            raw.append(row)
        # Move the file contents to out buffer list
        captureRate = int(raw[1][0])
        # Capture rate is in the first column of the second row of the CSV
        names = raw[2]
        # The joint names are located in the third row
        rows = raw[5:-1]
        # Read in data starting from the 6th row of the CSV

        for i, name in enumerate(names):
            if i%3 == 2:
                idx = name.rfind(':')
                name = name[idx+1:]
                jointNames.append(name)
            # The joint names begin at index 2 of the names list and appear at every third index after that.
            # The subject name and marker name are divided by a :
        # Create list of jointnames from label
        print("The available markers are ", jointNames)
        # Print the list of marker names to console so users know what they can choose to simulate.
        for row in rows:
            # print(row)
            timestamps.append(float(row[0]))
            # Set timestamp as first element of row
            newData = [float(element) for element in row[2:]]
            newData = np.split(np.array(newData), len(jointNames))
            data.append(newData)
            # Split data into x-y-z tuples to be processed by the marker capture object

    timestamps = np.array(timestamps)
    data = np.array(data)/1000.0
    # Data is in mm, so divide by 1000 to get m

    capture = MarkerCapture()
    capture.frameTimes = timestamps
    capture.frameCount = len(timestamps)
    capture.frameRate = captureRate
    capture.framePeriod = 1.0/captureRate
    # Create the marker capture object and set necessary attributes
    
    for i, name in enumerate(jointNames):
        marker = Marker3DOF(capture, name)
        for time, position in zip(capture.frameTimes, data[:,i,:]):
            marker.positionKeyFrames.add(time, position.reshape(3,1))
    # Set time and data attributes for each joint for the marker capture object
    return capture
    # Return the capture object that can then be splined and simulated