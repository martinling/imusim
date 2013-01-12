"""
Basestation hardware platform models.
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

from imusim.platforms.base import Platform
from imusim.platforms.timers import Timer, IdealTimer
from imusim.platforms.radios import Radio, IdealRadio

class IdealBasestation(Platform):
    """
    A simple basestation with ideal components.

    @ivar timer: The basestation L{Timer}.
    @ivar radio: The basestation L{Radio}
    """

    def __init__(self, simulation=None, trajectory=None):
        self.timer = IdealTimer(self)
        self.radio = IdealRadio(self)
        Platform.__init__(self, simulation, trajectory)

    @property
    def components(self):
        return [self.timer, self.radio]
