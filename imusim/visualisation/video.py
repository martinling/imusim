"""
Support for creating videos.
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

import numpy as np
from imusim.visualisation.rendering import AnimatedRenderer
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
import tempfile
import os
import subprocess
import shutil
import collections

def createVideo(filename, renderers, tMin, tMax, framerate=25,
        resolution=None):
    """
    Create a video of a 3D visualisation.

    The video is created using the current MayaVi figure, which should be setup
    beforehand to correctly position the camera and any other elements that
    should be rendered.

    @param filename: Name of the video file to create.
    @param renderers: List of L{AnimatedRenderer} objects to use.
    @param tMin: Time from which to begin rendering.
    @param tMax: Time at which to stop rendering.
    @param framerate: Framerate of the video in frames per second.
    @param resolution: Resolution of the video (width, height). If not
        specified, the resolution of the current figure is used.
    """
    if not isinstance(renderers, collections.Iterable):
        renderers = [renderers]

    figure = mlab.gcf()
    if resolution is not None:
        originalResolution = figure.scene.get_size()
        figure.scene.set_size(resolution)
    figure.scene.off_screen_rendering = True
    figure.scene.disable_render = True

    frameTimes = np.arange(tMin, tMax, 1.0/framerate)

    try:
        imageDir = tempfile.mkdtemp()
        for f,t in enumerate(frameTimes):
            for r in renderers:
                r.renderUpdate(t)
            framefile = os.path.join(imageDir, '%06d.png'%(f))
            mlab.savefig(framefile, size=resolution)
            assert os.path.exists(framefile)

        retval = subprocess.call(['ffmpeg', '-f','image2', '-r', str(framerate),
                '-i', '%s'%os.path.join(imageDir,'%06d.png'), '-sameq',
                filename, '-pass','2'])
    finally:
        shutil.rmtree(imageDir)
        figure.scene.disable_render = False
        figure.scene.off_screen_rendering = False
        if resolution is not None:
            figure.scene.set_size(originalResolution)
