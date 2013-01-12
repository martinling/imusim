"""
Sampling behaviours.
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

from imusim.platforms.adcs import ADC
from imusim.platforms.timers import Timer
from imusim.behaviours.timing import VirtualTimer
from imusim.platforms.sensors import Sensor

class PeriodicSampler(object):
    """
    Takes samples at regular intervals.
    """

    def __init__(self, adc, sensors, timer, period, handler):
        """
        Initialise sampler.

        The first sample will be taken after one period has elapsed.

        @param adc: L{ADC} to sample with.
        @param sensors: List of L{Sensor}s to sample.
        @param timer: L{Timer} or L{VirtualTimer}to use for timing.
        @param period: Interval between samples.
        @param handler: Function to call on each sample. The sampled
            values will be passed as arguments in the same order as the
            list of sensors passed to this constructor.
        """
        self._adc = adc
        self._sensors = sensors
        self._period = period
        self._handler = handler

        self._timer = timer

        timer.callback = self._timerHandler
        timer.start(period, repeat=True)

    def _timerHandler(self):
        self._adc.startSample(self._handler, *self._sensors)
