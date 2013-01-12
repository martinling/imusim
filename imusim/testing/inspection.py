"""
Utilities for code introspection.
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

import inspect

def getImplementations(module, superclass):
    """
    Get algorithm implementations from a module.

    @return: All concrete subclasses of the superclass found in the module.
    """
    predicate = lambda obj: inspect.isclass(obj) \
            and not inspect.isabstract(obj)
    return( imp for name,imp in inspect.getmembers(module,
            predicate) if issubclass(imp, superclass) )
