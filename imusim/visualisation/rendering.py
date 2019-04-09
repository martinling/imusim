"""
3D data rendering using Mayavi.
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

import time
import wx
import numpy as np
try:
    from mayavi import mlab
except ImportError:
    from enthought.mayavi import mlab
from abc import ABC, abstractmethod
try:
    from traits.api import HasTraits, Button, Range
    from traitsui.api import View, Group, Item
except ImportError:
    from enthought.traits.api import HasTraits, Button, Range
    from enthought.traits.ui.api import View, Group, Item
import imusim.maths.vectors as vectors


class InteractiveAnimation(HasTraits):
    """
    Animator that renders in the WxWidgets main loop, using idle events.
    """
    startButton = Button('Start Animation')
    stopButton = Button('Stop Animation')
    speed = Range(0.01, 1.0)

    traits_view = View(Group(Item('startButton'),
                             Item('stopButton'),
                             show_labels = False
                             ),
                            Item('_'),
                            Item(name = 'speed'),
                       title = 'Animation Controller',
                       buttons = ['OK'])

    def __init__(self, tMin, tMax, *renderers):
        """
        Initialise idle animator.

        @param tMin: Time to start animation.
        @param tMax: Time at which to wrap back to tMin.
        @param renderers: List of L{AnimatedRenderer} objects.
        """
        HasTraits.__init__(self)
        self._renderers = renderers
        self._tMin = tMin
        self._tMax = tMax
        self.speed = 1
        self.time = tMin
        self.playing = False
        self._scene = mlab.gcf().scene
        autoPositionCamera()
        for r in self._renderers:
            r.renderUpdate(self.time)
        self.edit_traits()

    def _startButton_fired(self):
        self.playing = True
        self.lastRenderTime = time.time()
        wx.GetApp().Bind(wx.EVT_IDLE, self._onIdleEvent)

    def start(self):
        """ Start animation. """
        self._startButton_fired()

    def _stopButton_fired(self):
        wx.GetApp().Bind(wx.EVT_IDLE, None)
        self.playing = False

    def stop(self):
        """ Stop animation. """
        self._stopButton_fired()

    def _onIdleEvent(self, event):
        if self.playing:
            self._scene.disable_render = True
            timeNow = time.time()
            t = self.time + (timeNow - self.lastRenderTime) * self.speed
            if t > self._tMax:
                self.time = self._tMin + t - self._tMax
            else:
                self.time = t
            self.time = t % self._tMax
            self.lastRenderTime = timeNow
            for r in self._renderers:
                r.renderUpdate(self.time)
            self._scene.disable_render = False
            event.RequestMore()

def autoPositionCamera():
    """
    Set camera position looking along the x-axis.
    """
    s = mlab.gcf().scene
    s.disable_render = True
    mlab.view(90,90,distance='auto',focalpoint='auto')
    mlab.roll(0)
    s.disable_render = False


class AnimatedRenderer(ABC):
    """
    Base class for animation renderers.
    """
    def __init__(self, *args, **kwargs):
        """
        Initialise renderer.
        """
        self._mlab_args = args
        self._mlab_kwargs = kwargs
        self._datasources = None

    @abstractmethod
    def renderSnapshot(self, t):
        """
        Render a snapshot a time `t`.

        @param t: Time(s) at which to render.

        @return: Mayavi sources for the rendered data.
        """

    @abstractmethod
    def renderUpdate(self, t):
        """
        Update the state of the Mayavi datasources at time t.
        """

    @property
    def datasources(self):
        """
        The Mayavi data sources that should be updated to perform animation.
        """
        if self._datasources is None:
            self._datasources = self.renderSnapshot(0)
        return self._datasources


class BodyModelRenderer(AnimatedRenderer):
    """
    Renders a body model using lines connecting points.

    Additional arguments are passed to the L{mayavi.mlab.plot3d}
    function to control the appearance of the rendered model.

    @param model: Root joint of the body model tree to render.
    """

    def __init__(self, model, *args, **kwargs):
        AnimatedRenderer.__init__(self, *args, **kwargs)
        self._endEffectors = [p for p in model.points
                if (not p.isJoint) or (not p.hasChildren)]

    def _generateDataPoints(self, t, endEffector):
        return np.hstack(j.position(t) for j in endEffector.ascendTree()
            if j.parent is None or vectors.norm(j.positionOffset) > 0)

    def renderSnapshot(self, t):
        sources = []
        for e in self._endEffectors:
            px, py, pz = self._generateDataPoints(t, e)
            sources.append(mlab.plot3d(px, py, pz, *self._mlab_args,
                **self._mlab_kwargs).mlab_source)
        return sources

    def renderUpdate(self, t):
        for e,s in zip(self._endEffectors,self.datasources):
            x,y,z = self._generateDataPoints(t, e)
            s.set(x=x, y=y, z=z)


class VelocityVectorRenderer(AnimatedRenderer):
    """
    Renders the velocity vectors of a trajectory.

    Additional arguments are passed to the L{mayavi.mlab.quiver3d}
    function to control the appearance of the rendered model.

    @param trajectory: The trajectory to render the velocity of.
    """

    def __init__(self, trajectory, *args, **kwargs):
        AnimatedRenderer.__init__(self, *args, **kwargs)
        self._trajectory = trajectory

    def _generateDataVector(self, t):
        px,py,pz = self._trajectory.position(t)
        vx,vy,vz = self._trajectory.velocity(t)
        return px,py,pz,vx,vy,vz

    def renderSnapshot(self, t):
        px,py,pz,vx,vy,vz = self._generateDataVector(t)
        return mlab.quiver3d(px,py,pz,vx,vy,vz, scale_factor=1,
                *self._mlab_args, **self._mlab_kwargs).mlab_source

    def renderUpdate(self, t):
        x,y,z,u,v,w = self._generateDataVector(t)
        self.datasources.set(x=x, y=y, z=z, u=u, v=v, w=w)


class ContactProbabilityRenderer(AnimatedRenderer):
    """
    Renders the contact probabilities of body model points.

    Intended to visualise results of the L{ContactTrackingKalmanFilter}.

    Additional arguments are passed to the L{mayavi.mlab.points3d}
    function to control the appearance of the rendered model.

    @param model: Body model to render contact probabilities for.
    """
    def __init__(self, model, *args, **kwargs):
        AnimatedRenderer.__init__(self, *args, **kwargs)
        self._points = [p for p in model.points
                if hasattr(p,'contactProbabilities')]

    def _generateData(self, t, point):
        px,py,pz = point.position(t)
        probability = point.contactProbabilities(t)
        return px,py,pz,probability

    def renderSnapshot(self, t):
        sources = []
        for j in self._points:
            px,py,pz,probability = self._generateData(t, j)
            sources.append(mlab.points3d(px,py,pz,probability,
                vmax=1, vmin=0, scale_mode='none', scale_factor=0.1,
                *self._mlab_args,
                **self._mlab_kwargs).mlab_source)

        return sources

    def renderUpdate(self, t):
        for s,j in zip(self.datasources, self._points):
            x,y,z,scalars = self._generateData(t, j)
            s.set(x=x, y=y, z=z, scalars=scalars)


def renderSolenoidCoil(solenoid, segments=1000, **kwargs):
    """
    Visualise L{SolenoidMagneticField} model coil in 3D.

    @param solenoid: L{SolenoidMagneticField} to render.
    @param segments: Number of coil segments to generate. Larger numbers
        result in smoother coils.
    @param kwargs: Additional keywords to pass to
        L{mayavi.mlab.plot3d}.
    """
    coils = 2*solenoid.b * solenoid.n

    z = np.linspace(-solenoid.b,solenoid.b,segments)
    theta = np.linspace(0,coils*2*np.pi,segments)
    x = solenoid.a * np.cos(theta)
    y = solenoid.a * np.sin(theta)

    x,y,z = solenoid.transform.apply(np.vstack((x,y,z)))

    if not 'tube_radius' in kwargs:kwargs['tube_radius']=0.001
    if not 'opacity' in kwargs:kwargs['opacity']=0.5
    mlab.plot3d(x,y,z,**kwargs)
