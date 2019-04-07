"""
Models of hardware timers.
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
from abc import abstractmethod
from imusim.platforms.base import Component
from imusim.simulation.base import Simulation
from imusim.utilities.documentation import prepend_method_doc
import SimPy.Simulation
import numpy as np


class Timer(Component):
    """
    Base class for simulated hardware timers.

    @ivar callback: Function to be called when timer fires.
    """
    def __init__(self, platform):
        self._process = None
        Component.__init__(self, platform)

    class _TimerProcess(SimPy.Simulation.Process):
        def __init__(self, timer):
            self.timer = timer
            SimPy.Simulation.Process.__init__(self,
                    sim=self.timer.platform.simulation.engine)
            self._lastTime = self.sim.now()
            self.sim.activate(self, self.eventLoop())

        def eventLoop(self):
            while True:
                yield SimPy.Simulation.hold, self, self.timer._period
                if self.timer._process is self:
                    if self.timer.callback is not None:
                        self.timer.callback()
                    if not self.timer._repeat:
                        self.timer._process = None
                        return
                else:
                    return

        def timeElapsed(self):
            return self.sim.now() - self._lastTime

    def _simulationChange(self):
        self._process = None

    def start(self, period, repeat=False):
        """
        Start the timer.

        @param period: Duration after which the timer should fire.
        @param repeat: If True, repeat at the given period, else fire once.
        """
        self._period = self._truePeriod(period)
        self._repeat = repeat
        self._process = self._TimerProcess(self)

    def clear(self):
        """
        Clear the timer.
        """
        self._process = None

    def timeElapsed(self):
        """
        Get the measured time elapsed since this timer was started.
        """
        return self._measuredPeriod(self._process.timeElapsed())

    @abstractmethod
    def _truePeriod(self, targetPeriod):
        """
        Convert an requested period to an actual elapsed period of the timer.
        """
        pass

    @abstractmethod
    def _measuredPeriod(self, _truePeriod):
        """
        Convert an actual elapsed period to a period reported by the timer.
        """
        pass


class IdealTimer(Timer):
    """
    An ideal timer with infinite resolution and no clock error.
    """

    def _truePeriod(self, targetPeriod):
        return targetPeriod

    def _measuredPeriod(self, _truePeriod):
        return _truePeriod


class ParametricTimer(Timer):
    """
    An timer with parametric resolution and frequency error.

    All period durations will be quantised by rounding up to the next multiple
    of the nominal period.

    The clock is modelled as having a fixed frequency offset which is randomly
    selected at initialisation time based on the freqError parameter.
    """

    @prepend_method_doc(Component)
    def __init__(self, platform, frequency, freqError, rng=None):
        """
        @param frequency: Clock frequency in Hz.
        @param freqError: Maximum frequency error in parts per million (ppm).
            The frequency error is assumed to be normally distributed with
            freqError corresponding to 3 standard deviations.
        @param rng: L{np.random.RandomState} from which to draw frequency
            error.
        """
        Timer.__init__(self, platform)

        if rng is None:
            rng = np.random.RandomState()

        self._idealFreq = frequency
        freqSigma = (freqError * 1E-6) / 3
        freqError = rng.normal(loc=1, scale=freqSigma)
        self._actualFreq = self._idealFreq * freqError

    def _truePeriod(self, targetPeriod):
        return np.ceil(targetPeriod * self._idealFreq) / self._actualFreq

    def _measuredPeriod(self, _truePeriod):
        return np.floor(_truePeriod * self._actualFreq) / self._idealFreq
