"""
Timing behaviours.
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
from abc import ABC, abstractmethod
from imusim.platforms.timers import Timer
from imusim.simulation.base import Simulation
from imusim.utilities.documentation import prepend_method_doc
import SimPy.Simulation
import numpy as np


class VirtualTimer(ABC):
    """
    A virtual timer multiplexed to a hardware timer.

    Implements the same interface as L{Timer}.

    @ivar callback: Function to be called when timer fires.
    """
    def __init__(self, mux, callback=None):
        """
        Create virtual timer.

        @param mux: The L{TimerMultiplexer} implementing this timer.
        """
        self._period = None
        self._mux = mux
        self.callback = callback
        self._mux._vtimers.append(self)

    @prepend_method_doc(Timer)
    def start(self, period, repeat=False):
        self._period = period
        self._repeat = repeat

        if self._mux._inHandler:
            # In timer event, handler is iterating through virtual timer list.
            self._done = True
            self._remaining = self._period
        else:
            self._done = False
            if self._mux._period is None:
                # Timer off, this will become the only active virtual timer.
                self._remaining = self._period
                self._mux._period = self._period
                self._mux._timer.start(self._period, self._mux._timerHandler)
            elif self._period < self._mux._period:
                # This timer will fire before any others.
                self._remaining = self._period
                elapsed = self._mux._timer.timeElapsed()
                for vtimer in self._mux._vtimers:
                    if vtimer is not self and vtimer._period is not None:
                        vtimer._remaining -= elapsed
                self._mux._period = self._period
                self._mux._timer.start(self._mux._period,
                        self._mux._timerHandler)
            elif self._period == self._mux._period:
                # This timer will fire at the same time as others.
                self._remaining = self._period
            else:
                # Other timers will fire before this one.
                self._remaining = self._period + self._mux._timer.timeElapsed()

    @prepend_method_doc(Timer)
    def clear(self):
        self._period = None

    @prepend_method_doc(Timer)
    def timeElapsed(self):
        return (self._period - self._remaining) + self.mux._timer.timeElapsed()


class TimerMultiplexer(ABC):
    """
    Multiplexer for running multiple virtual timers on one hardware timer.
    """
    def __init__(self, timer):
        """
        Create timer multiplexer.

        @param timer: L{Timer} to multiplex.
        """
        self._timer = timer
        self._vtimers = []
        self._period = None
        self._inHandler = False
        timer.clear()
        timer.callback = self._timerHandler

    def _timerHandler(self):
        """
        Timer event handler.
        """

        self._inHandler = True
        for vtimer in self._vtimers:
            if vtimer._period is not None and not vtimer._done:
                vtimer._remaining -= self._period
                if vtimer._remaining <= 0:
                    if vtimer._repeat:
                        vtimer._remaining = vtimer._period
                    else:
                        vtimer._period = None
                    if vtimer.callback is not None:
                        vtimer.callback()
        self._inHandler = False

        for vtimer in self._vtimers:
            vtimer._done = False

        waits = [vtimer._remaining for vtimer in self._vtimers
                if vtimer._period is not None]

        if len(waits) > 0:
            self._period = min(waits)
            self._timer.start(self._period, self._timerHandler)
        else:
            self._period = None
