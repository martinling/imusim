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
import functools


_rotationMatrices = dict(
    x = lambda rx: np.matrix((
        (1,0,0),
        (0,math.cos(rx),-math.sin(rx)),
        (0,math.sin(rx),math.cos(rx))),dtype=float),
    y = lambda ry: np.matrix((
        (math.cos(ry),0,math.sin(ry)),
        (0,1,0),
        (-math.sin(ry),0,math.cos(ry))),dtype=float),
    z = lambda rz: np.matrix((
        (math.cos(rz),-math.sin(rz),0),
        (math.sin(rz),math.cos(rz),0),
        (0,0,1)),dtype=float))

_EPS = 1e-12

_eulerFuncs = dict(
    xyz = lambda m: \
            (np.arctan2(-m[1,2], m[2,2]), np.arcsin(m[0,2]), np.arctan2(-m[0,1], m[0,0])) if abs(m[0,2]) < 1 - _EPS \
            else (np.arctan2(m[1,0], m[1,1]), np.pi/2, 0) if m[0,2] > 0 \
            else (-np.arctan2(m[1,0], m[1,1]), -np.pi/2, 0),
    xzy = lambda m: \
            (np.arctan2(m[2,1], m[1,1]), np.arcsin(-m[0,1]), np.arctan2(m[0,2], m[0,0])) if abs(m[0,1]) < 1 - _EPS \
            else (np.arctan2(-m[2,0], m[2,2]), -np.pi/2, 0) if m[0,1] > 0 \
            else (-np.arctan2(-m[2,0], m[2,2]), np.pi/2, 0),
    yxz = lambda m: \
            (np.arctan2(m[0,2], m[2,2]), np.arcsin(-m[1,2]), np.arctan2(m[1,0], m[1,1])) if abs(m[1,2]) < 1 - _EPS \
            else (np.arctan2(-m[0,1], m[0,0]), -np.pi/2, 0) if m[1,2] > 0 \
            else (-np.arctan2(-m[0,1], m[0,0]), np.pi/2, 0),
    yzx = lambda m: \
            (np.arctan2(-m[2,0], m[0,0]), np.arcsin(m[1,0]), np.arctan2(-m[1,2], m[1,1])) if abs(m[1,0]) < 1 - _EPS \
            else (np.arctan2(m[2,1], m[2,2]), np.pi/2, 0) if m[1,0] > 0 \
            else (-np.arctan2(m[2,1], m[2,2]), -np.pi/2, 0),
    zxy = lambda m: \
            (np.arctan2(-m[0,1], m[1,1]), np.arcsin(m[2,1]), np.arctan2(-m[2,0], m[2,2])) if abs(m[2,1]) < 1 - _EPS \
            else (np.arctan2(m[0,2], m[0,0]), np.pi/2, 0) if m[2,1] > 0 \
            else (-np.arctan2(m[0,2], m[0,0]), -np.pi/2, 0),
    zyx = lambda m: \
            (np.arctan2(m[1,0], m[0,0]), np.arcsin(-m[2,0]), np.arctan2(m[2,1], m[2,2])) if abs(m[2,0]) < 1 - _EPS \
            else (np.arctan2(-m[1,2], m[1,1]), -np.pi/2, 0) if m[2,0] > 0 \
            else (-np.arctan2(-m[1,2], m[1,1]), np.pi/2, 0),
    xyx = lambda m: \
            (np.arctan2(m[1,0], -m[2,0]), np.arccos(m[0,0]), np.arctan2(m[0,1], m[0,2])) if abs(m[0,0]) < 1 - _EPS \
            else (np.arctan2(-m[1,2], m[1,1]), 0, 0) if m[0,0] > 0 \
            else (-np.arctan2(-m[1,2], m[1,1]), np.pi, 0),
    xzx = lambda m: \
            (np.arctan2(m[2,0], m[1,0]), np.arccos(m[0,0]), np.arctan2(m[0,2], -m[0,1])) if abs(m[0,0]) < 1 - _EPS \
            else (np.arctan2(m[2,1], m[2,2]), 0, 0) if m[0,0] > 0 \
            else (-np.arctan2(m[2,1], m[2,2]), np.pi, 0),
    yxy = lambda m: \
            (np.arctan2(m[0,1], m[2,1]), np.arccos(m[1,1]), np.arctan2(m[1,0], -m[1,2])) if abs(m[1,1]) < 1 - _EPS \
            else (np.arctan2(m[0,2], m[0,0]), 0, 0) if m[1,1] > 0 \
            else (-np.arctan2(m[0,2], m[0,0]), np.pi, 0),
    yzy = lambda m: \
            (np.arctan2(m[2,1], -m[0,1]), np.arccos(m[1,1]), np.arctan2(m[1,2], m[1,0])) if abs(m[1,1]) < 1 - _EPS \
            else (np.arctan2(-m[2,0], m[2,2]), 0, 0) if m[1,1] > 0 \
            else (-np.arctan2(-m[2,0], m[2,2]), np.pi, 0),
    zxz = lambda m: \
            (np.arctan2(m[0,2], -m[1,2]), np.arccos(m[2,2]), np.arctan2(m[2,0], m[2,1])) if abs(m[2,2]) < 1 - _EPS \
            else (np.arctan2(-m[0,1], m[0,0]), 0, 0) if m[2,2] > 0 \
            else (-np.arctan2(-m[0,1], m[0,0]), np.pi, 0),
    zyz = lambda m: \
            (np.arctan2(m[1,2], m[0,2]), np.arccos(m[2,2]), np.arctan2(m[2,1], -m[2,0])) if abs(m[2,2]) < 1 - _EPS \
            else (np.arctan2(m[1,0], m[1,1]), 0, 0) if m[2,2] > 0 \
            else (-np.arctan2(m[1,0], m[1,1]), np.pi, 0),
    xy = lambda m: (np.arctan2(m[2,1], m[1,1]), np.arctan2(m[0,2], m[0,0])),
    xz = lambda m: (np.arctan2(-m[1,2], m[2,2]), np.arctan2(-m[0,1], m[0,0])),
    yx = lambda m: (np.arctan2(-m[2,0], m[0,0]), np.arctan2(-m[1,2], m[1,1])),
    yz = lambda m: (np.arctan2(m[0,2], m[2,2]), np.arctan2(m[1,0], m[1,1])),
    zx = lambda m: (np.arctan2(m[1,0], m[0,0]), np.arctan2(m[2,1], m[2,2])),
    zy = lambda m: (np.arctan2(-m[0,1], m[1,1]), np.arctan2(-m[2,0], m[2,2])),
    x = lambda m: (np.arctan2(m[2,1], m[2,2]),),
    y = lambda m: (np.arctan2(m[0,2], m[0,0]),),
    z = lambda m: (np.arctan2(m[1,0], m[1,1]),))

def matrixToEuler(m,order='zyx',inDegrees=True):
    """
    Convert a 3x3 rotation matrix to an Euler angle sequence.

    @param m: 3x3 L{np.matrix}, or equivalent, to convert.
    @param order: The order of the Euler angle sequence, e.g. 'zyx'
    @param inDegrees: True to return result in degrees, False for radians.

    @return: L{np.ndarray} of Euler angles in specified order.
    """

    order = order.lower()
    if order not in _eulerFuncs.keys():
        raise NotImplementedError("Order %s not implemented" % order)
    result = np.array(_eulerFuncs[order](m))

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
        angles = np.radians(angles)

    return functools.reduce(operator.mul,
            (_rotationMatrices[axis](angle) for axis,angle in
                zip(order.lower(), angles)))
