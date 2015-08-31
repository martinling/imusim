"""
Utilities for working with matrices.
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
import numpy as np
import math
import operator

_validEulerSequences = [
        'xyz', 'xzy',
        'yxz', 'yzx',
        'zxy', 'zyx',
        'xyx', 'xzx',
        'yxy', 'yzy',
        'zxz', 'zyz' ]

_rotationMatrices = dict(
    x = lambda rx: np.matrix((
        (1,0,0),
        (0,math.cos(rx),-math.sin(rx)),
        (0,math.sin(rx),math.cos(rx)))),
    y = lambda ry: np.matrix((
        (math.cos(ry),0,math.sin(ry)),
        (0,1,0),
        (-math.sin(ry),0,math.cos(ry)))),
    z = lambda rz: np.matrix((
        (math.cos(rz),-math.sin(rz),0),
        (math.sin(rz),math.cos(rz),0),
        (0,0,1))))

def matrixToEuler(m,order='zyx',inDegrees=True):
    """
    Convert a 3x3 rotation matrix to an Euler angle sequence.

    @param m: 3x3 L{np.matrix}, or equivalent, to convert.
    @param order: The order of the Euler angle sequence, 'zyx' or 'zxy'.
    @param inDegrees: True to return result in degrees, False for radians.

    @return: L{np.ndarray} of Euler angles in specified order.
    """
    EPS = 1e-6
    order = order.lower()
    assert order in _validEulerSequences, "Invalid Euler sequence '%s'" % order

    if order in 'zyx':
        sp = -m[2,0]
        if sp < (1-EPS):
            if sp > (-1+EPS):
                p = np.arcsin(sp)
                r = np.arctan2(m[2,1],m[2,2])
                y = np.arctan2(m[1,0],m[0,0])
            else:
                p = -np.pi/2.
                r = 0
                y = np.pi-np.arctan2(-m[0,1],m[0,2])
        else:
            p = np.pi/2.
            y = np.arctan2(-m[0,1],m[0,2])
            r = 0
        result = np.array((y,p,r))

    elif order in 'zxy':
        sx = m[2,1]
        if sx < (1-EPS):
            if sx > (-1+EPS):
                x = np.arcsin(sx)
                z = np.arctan2(-m[0,1],m[1,1])
                y = np.arctan2(-m[2,0],m[2,2])
            else:
                x = -np.pi/2
                y = 0
                z = -np.arctan2(m[0,2],m[0,0])
        else:
            x = np.pi/2
            y = 0
            z = np.arctan2(m[0,2],m[0,0])
        result = np.array((z,x,y))

    else:
        raise NotImplementedError, "Unimplemented rotation order:%r"%order

    if inDegrees:
        return np.degrees(result)
    else:
        return result

def matrixFromEuler(angles, order, inDegrees=True):
    """
    Generate a rotation matrix from an Euler angle sequence.

    @param angles: Sequence of Euler rotation angles.
    @param order: Sequence of rotation axes. Rotations are applied sequentially
        from left to right, i.e. the string 'zyx' would result in rotation about
        the Z axis, then the new Y axis, and finally about the new X axis.
    @param inDegrees: Whether the angles are in degrees (`True`) or radians
        (`False`)
    @return: 3x3 rotation matrix corresponding to the Euler angle sequence.
    """

    assert len(angles) == len(order)

    if inDegrees:
        angles = (np.radians(angle) for angle in angles)

    return reduce(operator.mul,
            (_rotationMatrices[axis.lower()](angle) for axis,angle in
                zip(order, angles)))



