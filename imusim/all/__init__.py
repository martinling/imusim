"""
Package which imports all commonly used symbols from IMUSim.
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

import pkgutil
import os
import inspect

__all__ = []
path = os.path.split(os.path.dirname(pkgutil.get_loader('imusim').path))[0]
for loader, modname, ispkg in pkgutil.walk_packages([path]):
    if modname in ['imusim.visualisation.rendering', 'imusim.visualisation.video']:
        continue
    if modname.startswith('imusim') \
            and not modname.startswith('imusim.tests') \
            and not modname.startswith('imusim.all'):
        exec("import %s" % modname)
        exec("module = %s" % modname)
        symbols = filter(lambda o: not inspect.ismodule(o),
                module.__all__ if hasattr(module,'__all__') else dir(module))
        symbols = filter(lambda s: not s.startswith('_'), symbols)
        exec("from %s import *" % modname)
        __all__ += symbols
