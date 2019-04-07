"""
Base class for simulated hardware platforms.
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
from imusim.simulation.base import Simulation

class Platform(ABC):
    """
    Base class for simulated hardware platforms.
    """
    def __init__(self, simulation=None, trajectory=None):
        """
        Initialise simulated hardware platform.

        @param simulation: The L{Simulation} this platform will be used in.
        @param trajectory: The trajectory to be followed by this platform.
        """
        self._simulation = simulation
        self._trajectory = trajectory
        self._simulationChange()
        self._trajectoryChange()

    @property
    @abstractmethod
    def components(self):
        """
        List of the components attached to this platform.
        """
        pass

    @property
    def simulation(self):
        """
        The L{Simulation} this platform is part of.
        """
        return self._simulation

    @simulation.setter
    def simulation(self, simulation):
        self._simulation = simulation
        self._simulationChange()

    @property
    def trajectory(self):
        """
        The trajectory followed by this platform.
        """
        return self._trajectory

    @trajectory.setter
    def trajectory(self, trajectory):
        self._trajectory = trajectory
        self._trajectoryChange()

    def _simulationChange(self):
        """
        Called when the platform is assigned a new simulation.
        """
        for component in self.components:
            component._simulationChange()

    def _trajectoryChange(self):
        """
        Called when the platform is assigned a new trajectory.
        """
        for component in self.components:
            component._trajectoryChange()


class Component(ABC):
    """
    Base class for simulated hardware components.
    """
    def __init__(self, platform):
        """
        Initialise simulated component.

        @param platform: The L{Platform} this component will be attached to.
        """
        self.platform = platform

    def _simulationChange(self):
        """
        Called when the platform is assigned a new simulation.
        """
        pass

    def _trajectoryChange(self):
        """
        Called when the platform is assigned a new trajectory.
        """
        pass
