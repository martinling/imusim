"""
Simulation engine for IMUSim.
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

import SimPy.Simulation
import numpy as np
import time

from imusim.environment.base import Environment
from imusim.maths.quaternions import Quaternion

class Simulation(object):
    """
    Main IMUSim simulation environment.

    @ivar rng: L{np.random.RandomState} instance for random number generation.
    @ivar environment: L{Environment} in use for simulation.
    @ivar engine: L{SimPy.Simulation} instance driving simulation.
    """

    def __init__(self, seed=None, environment=None,
            engine=SimPy.Simulation.Simulation):
        """
        Initialise simulation.

        @param seed: Seed value for L{np.random.RandomState}.
        @param environment: L{Environment} to simulate in.
        @param engine: L{SimPy.Simulation} subclass to drive simulation.
        """
        self.environment = Environment() if environment is None else environment
        self.engine = engine()
        self.engine.initialize()
        self.rng = np.random.RandomState(seed)

    @property
    def time(self):
        """ The current true time within the simulation (float). """
        assert hasattr(self, '_startTime'), 'Simulation time has not been set.'
        return self.engine.now() + self._startTime

    @time.setter
    def time(self, time):
        assert not hasattr(self, '_startTime'), 'Simulation time already set.'
        self._startTime = time

    class ProgressMonitor(SimPy.Simulation.Process):
        """
        Simple process to print status and remaining simulation time at 5%
        progress intervals.
        """
        def __init__(self, sim, duration):
            SimPy.Simulation.Process.__init__(self, sim=sim)
            self.startTime = sim.now()
            self.duration = duration
            sim.activate(self, self.printProgress())

        def printProgress(self):
            for i in range(0,20):
                start = time.time()
                yield SimPy.Simulation.hold, self, self.duration / 20.0
                now = self.sim.now() - self.startTime
                end = time.time()

                remaining = (20-i) * (end-start)
                print ("Simulated %.1fs of %.1fs (%3.0f%%). Estimated time remaining %.1fs" % (
                    now, self.duration, (i + 1) * 5, remaining))

    def run(self, endTime, printProgress=True):
        """
        Advance the simulation up to the given time.

        @param endTime: Simulation time at which to stop (float).
        @param printProgress: Whether to print progress information (bool).
        """

        if printProgress:
            print("Simulating...")
            startTime = self.time
            duration = endTime - startTime
            progressMonitor = self.ProgressMonitor(self.engine, duration)

        startWallTime = time.time()
        self.engine.simulate(endTime - self._startTime)
        endWallTime = time.time()

        if printProgress:
            print("Simulation complete.")
            print("Simulated %.1f seconds in %.1f seconds." % (
                    endTime - startTime, endWallTime - startWallTime))

    def subrng(self):
        """
        Obtain a new random number generator seeded from the main RNG.

        @return: A L{np.random.RandomState} instance.
        """
        return np.random.RandomState(seed=int(self.rng.uniform(0,2**31)))
